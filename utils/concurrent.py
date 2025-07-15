"""
Concurrent processing utilities for bioinformatics operations.

This module provides rate-limited executors, batch processors, and progress
tracking for parallel bioinformatics computations.
"""

import time
import threading
from typing import Any, Callable, Dict, List, Optional, Union, Iterator
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed, Future
from dataclasses import dataclass
from enum import Enum
import multiprocessing as mp
import queue
import logging
from pathlib import Path
import psutil

from ..core.exceptions import ProcessingError
from ..core.config import get_config


class ExecutorType(Enum):
    """Types of parallel executors."""
    THREAD = "thread"
    PROCESS = "process"


@dataclass
class TaskResult:
    """Result of a parallel task."""
    task_id: str
    result: Any
    success: bool
    error: Optional[Exception] = None
    execution_time: float = 0.0
    memory_usage: Optional[float] = None


class ProgressTracker:
    """
    Thread-safe progress tracker for long-running operations.
    
    Provides real-time progress updates and performance statistics.
    """
    
    def __init__(self, total_tasks: int, description: str = "Processing"):
        self.total_tasks = total_tasks
        self.description = description
        self.completed_tasks = 0
        self.failed_tasks = 0
        self.start_time = time.time()
        self.lock = threading.Lock()
        self.logger = logging.getLogger(__name__)
        
        # Performance tracking
        self.task_times: List[float] = []
        self.memory_usage: List[float] = []
    
    def update(self, increment: int = 1, task_time: Optional[float] = None, memory_mb: Optional[float] = None):
        """Update progress with optional performance data."""
        with self.lock:
            self.completed_tasks += increment
            
            if task_time is not None:
                self.task_times.append(task_time)
            
            if memory_mb is not None:
                self.memory_usage.append(memory_mb)
            
            # Log progress at intervals
            if self.completed_tasks % max(1, self.total_tasks // 10) == 0:
                self._log_progress()
    
    def mark_failed(self, increment: int = 1):
        """Mark tasks as failed."""
        with self.lock:
            self.failed_tasks += increment
            self.completed_tasks += increment
            
            if self.completed_tasks % max(1, self.total_tasks // 10) == 0:
                self._log_progress()
    
    def _log_progress(self):
        """Log current progress."""
        percentage = (self.completed_tasks / self.total_tasks) * 100
        elapsed = time.time() - self.start_time
        
        if self.completed_tasks > 0:
            avg_time = elapsed / self.completed_tasks
            eta = avg_time * (self.total_tasks - self.completed_tasks)
            
            message = (
                f"{self.description}: {self.completed_tasks}/{self.total_tasks} "
                f"({percentage:.1f}%) - ETA: {eta:.1f}s"
            )
            
            if self.failed_tasks > 0:
                message += f" - Failed: {self.failed_tasks}"
            
            self.logger.info(message)
    
    def get_stats(self) -> Dict[str, Any]:
        """Get comprehensive progress statistics."""
        with self.lock:
            elapsed = time.time() - self.start_time
            
            stats = {
                'total_tasks': self.total_tasks,
                'completed_tasks': self.completed_tasks,
                'failed_tasks': self.failed_tasks,
                'success_rate': (self.completed_tasks - self.failed_tasks) / max(1, self.completed_tasks),
                'elapsed_time': elapsed,
                'completion_percentage': (self.completed_tasks / self.total_tasks) * 100
            }
            
            if self.task_times:
                stats.update({
                    'avg_task_time': sum(self.task_times) / len(self.task_times),
                    'min_task_time': min(self.task_times),
                    'max_task_time': max(self.task_times),
                    'total_processing_time': sum(self.task_times)
                })
            
            if self.memory_usage:
                stats.update({
                    'avg_memory_mb': sum(self.memory_usage) / len(self.memory_usage),
                    'max_memory_mb': max(self.memory_usage),
                    'min_memory_mb': min(self.memory_usage)
                })
            
            return stats
    
    def is_complete(self) -> bool:
        """Check if all tasks are complete."""
        with self.lock:
            return self.completed_tasks >= self.total_tasks


class RateLimitedExecutor:
    """
    Rate-limited executor for API calls and external services.
    
    Provides controlled execution with rate limiting, retry logic,
    and adaptive delays.
    """
    
    def __init__(
        self,
        rate_limit: float = 0.2,  # Seconds between requests
        max_workers: int = 4,
        executor_type: ExecutorType = ExecutorType.THREAD,
        retry_attempts: int = 3,
        backoff_factor: float = 2.0
    ):
        self.rate_limit = rate_limit
        self.max_workers = max_workers
        self.executor_type = executor_type
        self.retry_attempts = retry_attempts
        self.backoff_factor = backoff_factor
        
        self.last_request_time = 0.0
        self.request_lock = threading.Lock()
        self.logger = logging.getLogger(__name__)
        
        # Create appropriate executor
        if executor_type == ExecutorType.THREAD:
            self.executor = ThreadPoolExecutor(max_workers=max_workers)
        else:
            self.executor = ProcessPoolExecutor(max_workers=max_workers)
    
    def _apply_rate_limit(self):
        """Apply rate limiting between requests."""
        with self.request_lock:
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            
            if time_since_last < self.rate_limit:
                sleep_time = self.rate_limit - time_since_last
                time.sleep(sleep_time)
            
            self.last_request_time = time.time()
    
    def _execute_with_retry(self, func: Callable, *args, **kwargs) -> Any:
        """Execute function with retry logic."""
        last_exception = None
        
        for attempt in range(self.retry_attempts):
            try:
                self._apply_rate_limit()
                return func(*args, **kwargs)
                
            except Exception as e:
                last_exception = e
                if attempt < self.retry_attempts - 1:
                    delay = self.rate_limit * (self.backoff_factor ** attempt)
                    self.logger.warning(
                        f"Attempt {attempt + 1} failed, retrying in {delay:.2f}s: {e}"
                    )
                    time.sleep(delay)
                else:
                    self.logger.error(f"All {self.retry_attempts} attempts failed")
        
        raise last_exception
    
    def submit(self, func: Callable, *args, **kwargs) -> Future:
        """Submit function for rate-limited execution."""
        return self.executor.submit(self._execute_with_retry, func, *args, **kwargs)
    
    def map(self, func: Callable, iterable: Iterator, **kwargs) -> List[Any]:
        """Map function over iterable with rate limiting."""
        futures = [self.submit(func, item, **kwargs) for item in iterable]
        results = []
        
        for future in as_completed(futures):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                self.logger.error(f"Task failed: {e}")
                results.append(None)
        
        return results
    
    def shutdown(self, wait: bool = True):
        """Shutdown the executor."""
        self.executor.shutdown(wait=wait)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.shutdown()


class BatchProcessor:
    """
    High-performance batch processor for bioinformatics workflows.
    
    Provides memory-aware batch processing with progress tracking,
    error handling, and performance optimization.
    """
    
    def __init__(
        self,
        batch_size: int = 100,
        max_workers: int = None,
        executor_type: ExecutorType = ExecutorType.THREAD,
        memory_limit_gb: Optional[float] = None,
        progress_callback: Optional[Callable] = None
    ):
        self.batch_size = batch_size
        self.max_workers = max_workers or min(32, (mp.cpu_count() or 1) + 4)
        self.executor_type = executor_type
        self.memory_limit_gb = memory_limit_gb
        self.progress_callback = progress_callback
        self.logger = logging.getLogger(__name__)
        
        # Performance monitoring
        self.total_processed = 0
        self.total_failed = 0
        self.start_time = None
    
    def create_batches(self, items: List[Any], batch_size: Optional[int] = None) -> List[List[Any]]:
        """Split items into batches for processing."""
        batch_size = batch_size or self.batch_size
        batches = []
        
        for i in range(0, len(items), batch_size):
            batch = items[i:i + batch_size]
            batches.append(batch)
        
        return batches
    
    def _monitor_memory(self) -> float:
        """Monitor current memory usage in GB."""
        process = psutil.Process()
        memory_mb = process.memory_info().rss / 1024 / 1024
        return memory_mb / 1024
    
    def _process_batch_with_monitoring(
        self,
        batch_func: Callable,
        batch: List[Any],
        batch_id: str,
        **kwargs
    ) -> TaskResult:
        """Process batch with performance monitoring."""
        start_time = time.time()
        start_memory = self._monitor_memory()
        
        try:
            # Check memory limit
            if self.memory_limit_gb and start_memory > self.memory_limit_gb:
                raise ProcessingError(
                    f"Memory limit exceeded: {start_memory:.2f}GB > {self.memory_limit_gb}GB"
                )
            
            # Process batch
            result = batch_func(batch, **kwargs)
            
            # Calculate performance metrics
            end_time = time.time()
            end_memory = self._monitor_memory()
            execution_time = end_time - start_time
            memory_usage = max(end_memory - start_memory, 0)
            
            return TaskResult(
                task_id=batch_id,
                result=result,
                success=True,
                execution_time=execution_time,
                memory_usage=memory_usage
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            self.logger.error(f"Batch {batch_id} failed after {execution_time:.2f}s: {e}")
            
            return TaskResult(
                task_id=batch_id,
                result=None,
                success=False,
                error=e,
                execution_time=execution_time
            )
    
    def process_parallel(
        self,
        items: List[Any],
        batch_func: Callable,
        description: str = "Processing batches",
        **kwargs
    ) -> List[TaskResult]:
        """
        Process items in parallel batches with comprehensive monitoring.
        
        Args:
            items: List of items to process
            batch_func: Function to process each batch
            description: Progress description
            **kwargs: Additional arguments for batch_func
            
        Returns:
            List of TaskResult objects
        """
        if not items:
            return []
        
        self.start_time = time.time()
        batches = self.create_batches(items)
        
        # Initialize progress tracking
        tracker = ProgressTracker(len(batches), description)
        
        self.logger.info(
            f"Processing {len(items)} items in {len(batches)} batches "
            f"with {self.max_workers} workers"
        )
        
        # Choose executor type
        executor_class = ThreadPoolExecutor if self.executor_type == ExecutorType.THREAD else ProcessPoolExecutor
        
        results = []
        
        with executor_class(max_workers=self.max_workers) as executor:
            # Submit all batches
            future_to_batch = {}
            for i, batch in enumerate(batches):
                batch_id = f"batch_{i:04d}"
                future = executor.submit(
                    self._process_batch_with_monitoring,
                    batch_func,
                    batch,
                    batch_id,
                    **kwargs
                )
                future_to_batch[future] = (batch_id, batch)
            
            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch_id, batch = future_to_batch[future]
                
                try:
                    result = future.result()
                    results.append(result)
                    
                    if result.success:
                        self.total_processed += len(batch)
                        tracker.update(
                            1,
                            task_time=result.execution_time,
                            memory_mb=result.memory_usage * 1024 if result.memory_usage else None
                        )
                    else:
                        self.total_failed += len(batch)
                        tracker.mark_failed(1)
                    
                    # Call progress callback if provided
                    if self.progress_callback:
                        self.progress_callback(tracker.get_stats())
                        
                except Exception as e:
                    self.logger.error(f"Future failed for {batch_id}: {e}")
                    self.total_failed += len(batch)
                    tracker.mark_failed(1)
                    
                    # Add error result
                    error_result = TaskResult(
                        task_id=batch_id,
                        result=None,
                        success=False,
                        error=e
                    )
                    results.append(error_result)
        
        # Log final statistics
        total_time = time.time() - self.start_time
        success_rate = (self.total_processed / (self.total_processed + self.total_failed)) * 100
        
        self.logger.info(
            f"Batch processing complete: {self.total_processed} processed, "
            f"{self.total_failed} failed ({success_rate:.1f}% success) "
            f"in {total_time:.2f}s"
        )
        
        return results
    
    def process_sequential_batches(
        self,
        items: List[Any],
        batch_func: Callable,
        description: str = "Processing batches",
        **kwargs
    ) -> List[TaskResult]:
        """
        Process batches sequentially (useful for memory-intensive operations).
        
        Args:
            items: List of items to process
            batch_func: Function to process each batch
            description: Progress description
            **kwargs: Additional arguments for batch_func
            
        Returns:
            List of TaskResult objects
        """
        if not items:
            return []
        
        self.start_time = time.time()
        batches = self.create_batches(items)
        tracker = ProgressTracker(len(batches), description)
        
        results = []
        
        for i, batch in enumerate(batches):
            batch_id = f"batch_{i:04d}"
            result = self._process_batch_with_monitoring(batch_func, batch, batch_id, **kwargs)
            results.append(result)
            
            if result.success:
                self.total_processed += len(batch)
                tracker.update(
                    1,
                    task_time=result.execution_time,
                    memory_mb=result.memory_usage * 1024 if result.memory_usage else None
                )
            else:
                self.total_failed += len(batch)
                tracker.mark_failed(1)
            
            # Call progress callback if provided
            if self.progress_callback:
                self.progress_callback(tracker.get_stats())
        
        return results


class MemoryAwareBatchProcessor(BatchProcessor):
    """
    Memory-aware batch processor that adapts batch size based on memory usage.
    
    Automatically adjusts batch sizes to stay within memory limits and
    provides memory usage optimization.
    """
    
    def __init__(self, target_memory_gb: float = 4.0, **kwargs):
        super().__init__(**kwargs)
        self.target_memory_gb = target_memory_gb
        self.adaptive_batch_size = self.batch_size
        self.memory_samples = []
    
    def _adapt_batch_size(self, memory_usage_gb: float, processing_time: float):
        """Adapt batch size based on memory usage and processing time."""
        # Calculate memory efficiency (items per GB per second)
        if memory_usage_gb > 0 and processing_time > 0:
            efficiency = self.adaptive_batch_size / (memory_usage_gb * processing_time)
            
            # Adjust batch size to target memory usage
            if memory_usage_gb > self.target_memory_gb * 1.2:
                # Reduce batch size if using too much memory
                self.adaptive_batch_size = max(1, int(self.adaptive_batch_size * 0.8))
                self.logger.info(f"Reduced batch size to {self.adaptive_batch_size} (high memory usage)")
            elif memory_usage_gb < self.target_memory_gb * 0.5 and processing_time < 10:
                # Increase batch size if using little memory and processing quickly
                self.adaptive_batch_size = min(self.batch_size * 2, int(self.adaptive_batch_size * 1.2))
                self.logger.info(f"Increased batch size to {self.adaptive_batch_size} (low memory usage)")
    
    def create_batches(self, items: List[Any], batch_size: Optional[int] = None) -> List[List[Any]]:
        """Create batches with adaptive sizing."""
        batch_size = batch_size or self.adaptive_batch_size
        return super().create_batches(items, batch_size)
    
    def process_parallel(self, items: List[Any], batch_func: Callable, **kwargs) -> List[TaskResult]:
        """Process with adaptive memory management."""
        results = super().process_parallel(items, batch_func, **kwargs)
        
        # Analyze memory usage for future optimization
        memory_usages = [r.memory_usage for r in results if r.memory_usage]
        processing_times = [r.execution_time for r in results if r.success]
        
        if memory_usages and processing_times:
            avg_memory = sum(memory_usages) / len(memory_usages)
            avg_time = sum(processing_times) / len(processing_times)
            self._adapt_batch_size(avg_memory, avg_time)
        
        return results