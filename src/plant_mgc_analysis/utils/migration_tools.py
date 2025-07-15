"""
Migration utilities for Plant MGC Analysis Pipeline.

This module provides tools to migrate legacy scripts to the new unified architecture,
eliminating duplications and converting to OOP patterns.
"""

import re
import ast
import shutil
from typing import Dict, List, Optional, Tuple, Set
from pathlib import Path
from dataclasses import dataclass, field

from loguru import logger

from ..config.settings import get_settings
from ..utils.file_operations import file_manager


@dataclass
class ScriptAnalysis:
    """Analysis of a legacy script."""
    
    script_path: Path
    imports: List[str] = field(default_factory=list)
    functions: List[str] = field(default_factory=list)
    classes: List[str] = field(default_factory=list)
    hardcoded_paths: List[str] = field(default_factory=list)
    duplicated_patterns: List[str] = field(default_factory=list)
    complexity_score: int = 0
    migration_priority: str = "medium"
    suggested_module: str = "utils"
    
    def __post_init__(self):
        """Calculate complexity score."""
        self.complexity_score = (
            len(self.functions) * 2 +
            len(self.classes) * 5 +
            len(self.hardcoded_paths) * 3 +
            len(self.duplicated_patterns) * 4
        )
        
        # Determine migration priority
        if self.complexity_score > 50:
            self.migration_priority = "high"
        elif self.complexity_score > 20:
            self.migration_priority = "medium"
        else:
            self.migration_priority = "low"


@dataclass
class MigrationPlan:
    """Migration plan for legacy scripts."""
    
    total_scripts: int
    analyzed_scripts: List[ScriptAnalysis] = field(default_factory=list)
    high_priority: List[ScriptAnalysis] = field(default_factory=list)
    medium_priority: List[ScriptAnalysis] = field(default_factory=list)
    low_priority: List[ScriptAnalysis] = field(default_factory=list)
    duplications: Dict[str, List[Path]] = field(default_factory=dict)
    
    def __post_init__(self):
        """Organize scripts by priority."""
        self.high_priority = [s for s in self.analyzed_scripts if s.migration_priority == "high"]
        self.medium_priority = [s for s in self.analyzed_scripts if s.migration_priority == "medium"]
        self.low_priority = [s for s in self.analyzed_scripts if s.migration_priority == "low"]
    
    def get_summary(self) -> Dict[str, any]:
        """Get migration plan summary."""
        return {
            "total_scripts": self.total_scripts,
            "analyzed": len(self.analyzed_scripts),
            "high_priority": len(self.high_priority),
            "medium_priority": len(self.medium_priority),
            "low_priority": len(self.low_priority),
            "total_duplications": len(self.duplications),
            "estimated_effort_hours": self._estimate_effort(),
        }
    
    def _estimate_effort(self) -> int:
        """Estimate migration effort in hours."""
        return (
            len(self.high_priority) * 4 +
            len(self.medium_priority) * 2 +
            len(self.low_priority) * 1
        )


class LegacyScriptAnalyzer:
    """Analyzer for legacy scripts."""
    
    def __init__(self):
        """Initialize analyzer."""
        self.settings = get_settings()
        self.logger = logger
        
        # Common patterns to detect
        self.hardcoded_path_patterns = [
            r'/groups/itay_mayrose/alongonda/[^\'"\s]+',
            r'/tmp/[^\'"\s]+',
            r'/var/[^\'"\s]+',
            r'~/[^\'"\s]+',
        ]
        
        self.duplication_patterns = [
            r'pd\.read_csv\(',
            r'SeqIO\.parse\(',
            r'os\.path\.join\(',
            r'subprocess\.run\(',
            r'for record in SeqIO\.parse',
            r'with open\([^)]+\) as',
            r'csv\.DictReader\(',
            r'json\.load\(',
            r'requests\.get\(',
        ]
        
        # Function patterns that indicate specific functionality
        self.function_patterns = {
            'blast': [r'blast', r'makeblastdb', r'blastp', r'blastn'],
            'fasta': [r'fasta', r'SeqIO', r'sequence'],
            'csv': [r'csv', r'pandas', r'DataFrame'],
            'kegg': [r'kegg', r'pathway', r'rest\.kegg'],
            'analysis': [r'analyze', r'calculate', r'process'],
            'file_ops': [r'file', r'path', r'directory'],
        }
    
    def analyze_script(self, script_path: Path) -> ScriptAnalysis:
        """Analyze a single legacy script."""
        try:
            with open(script_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Parse AST
            try:
                tree = ast.parse(content)
            except SyntaxError as e:
                self.logger.warning(f"Syntax error in {script_path}: {e}")
                return ScriptAnalysis(script_path=script_path)
            
            # Extract information
            imports = self._extract_imports(tree)
            functions = self._extract_functions(tree)
            classes = self._extract_classes(tree)
            hardcoded_paths = self._find_hardcoded_paths(content)
            duplicated_patterns = self._find_duplicated_patterns(content)
            suggested_module = self._suggest_module(script_path, functions, imports)
            
            analysis = ScriptAnalysis(
                script_path=script_path,
                imports=imports,
                functions=functions,
                classes=classes,
                hardcoded_paths=hardcoded_paths,
                duplicated_patterns=duplicated_patterns,
                suggested_module=suggested_module,
            )
            
            self.logger.debug(f"Analyzed script: {script_path}")
            return analysis
            
        except Exception as e:
            self.logger.error(f"Error analyzing script {script_path}: {e}")
            return ScriptAnalysis(script_path=script_path)
    
    def _extract_imports(self, tree: ast.AST) -> List[str]:
        """Extract import statements."""
        imports = []
        
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imports.append(alias.name)
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                for alias in node.names:
                    imports.append(f"{module}.{alias.name}")
        
        return imports
    
    def _extract_functions(self, tree: ast.AST) -> List[str]:
        """Extract function definitions."""
        functions = []
        
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                functions.append(node.name)
        
        return functions
    
    def _extract_classes(self, tree: ast.AST) -> List[str]:
        """Extract class definitions."""
        classes = []
        
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                classes.append(node.name)
        
        return classes
    
    def _find_hardcoded_paths(self, content: str) -> List[str]:
        """Find hardcoded paths in content."""
        paths = []
        
        for pattern in self.hardcoded_path_patterns:
            matches = re.findall(pattern, content)
            paths.extend(matches)
        
        return list(set(paths))
    
    def _find_duplicated_patterns(self, content: str) -> List[str]:
        """Find duplicated code patterns."""
        patterns = []
        
        for pattern in self.duplication_patterns:
            if re.search(pattern, content):
                patterns.append(pattern)
        
        return patterns
    
    def _suggest_module(self, script_path: Path, functions: List[str], imports: List[str]) -> str:
        """Suggest appropriate module for script."""
        script_name = script_path.stem.lower()
        script_content = script_name + " " + " ".join(functions) + " " + " ".join(imports)
        
        # Check function patterns
        for module, patterns in self.function_patterns.items():
            if any(re.search(pattern, script_content, re.IGNORECASE) for pattern in patterns):
                return module
        
        # Check by directory
        if "blast" in str(script_path).lower():
            return "genomics"
        elif "kegg" in str(script_path).lower():
            return "metabolic"
        elif "analysis" in str(script_path).lower():
            return "analysis"
        else:
            return "utils"
    
    def analyze_directory(self, directory: Path) -> MigrationPlan:
        """Analyze all scripts in directory."""
        python_files = list(directory.glob("**/*.py"))
        
        # Filter out __init__.py and already migrated files
        python_files = [
            f for f in python_files
            if f.name != "__init__.py" and not f.name.endswith("_oop.py")
        ]
        
        self.logger.info(f"Analyzing {len(python_files)} Python files")
        
        analyses = []
        for script_path in python_files:
            analysis = self.analyze_script(script_path)
            analyses.append(analysis)
        
        # Find duplications
        duplications = self._find_duplications(analyses)
        
        plan = MigrationPlan(
            total_scripts=len(python_files),
            analyzed_scripts=analyses,
            duplications=duplications,
        )
        
        self.logger.info(f"Migration plan created: {plan.get_summary()}")
        return plan
    
    def _find_duplications(self, analyses: List[ScriptAnalysis]) -> Dict[str, List[Path]]:
        """Find duplicated functionality across scripts."""
        duplications = {}
        
        # Group by function names
        function_groups = {}
        for analysis in analyses:
            for func in analysis.functions:
                if func not in function_groups:
                    function_groups[func] = []
                function_groups[func].append(analysis.script_path)
        
        # Find duplicates
        for func, paths in function_groups.items():
            if len(paths) > 1:
                duplications[func] = paths
        
        return duplications


class ScriptMigrator:
    """Tool to migrate legacy scripts to new architecture."""
    
    def __init__(self):
        """Initialize migrator."""
        self.settings = get_settings()
        self.logger = logger
        self.analyzer = LegacyScriptAnalyzer()
    
    def migrate_script(
        self,
        script_path: Path,
        target_module: str,
        dry_run: bool = True
    ) -> Dict[str, any]:
        """Migrate a single script to new architecture."""
        # Analyze script
        analysis = self.analyzer.analyze_script(script_path)
        
        # Read original content
        with open(script_path, 'r', encoding='utf-8') as f:
            original_content = f.read()
        
        # Transform content
        transformed_content = self._transform_script_content(
            original_content, analysis, target_module
        )
        
        # Determine target path
        target_path = self._get_target_path(script_path, target_module)
        
        migration_result = {
            "original_path": str(script_path),
            "target_path": str(target_path),
            "analysis": analysis,
            "transformed": transformed_content is not None,
            "dry_run": dry_run,
        }
        
        if not dry_run and transformed_content:
            # Write transformed script
            target_path.parent.mkdir(parents=True, exist_ok=True)
            with open(target_path, 'w', encoding='utf-8') as f:
                f.write(transformed_content)
            
            self.logger.info(f"Migrated script: {script_path} -> {target_path}")
        
        return migration_result
    
    def _transform_script_content(
        self,
        content: str,
        analysis: ScriptAnalysis,
        target_module: str
    ) -> Optional[str]:
        """Transform script content to new architecture."""
        try:
            # Replace hardcoded paths
            content = self._replace_hardcoded_paths(content)
            
            # Add imports for new architecture
            content = self._add_new_imports(content, target_module)
            
            # Transform to OOP pattern
            content = self._transform_to_oop(content, analysis)
            
            # Add proper error handling
            content = self._add_error_handling(content)
            
            # Add logging
            content = self._add_logging(content)
            
            return content
            
        except Exception as e:
            self.logger.error(f"Error transforming script content: {e}")
            return None
    
    def _replace_hardcoded_paths(self, content: str) -> str:
        """Replace hardcoded paths with configuration."""
        path_mapping = self.settings.migrate_legacy_paths()
        
        for old_path, new_path in path_mapping.items():
            content = content.replace(f'"{old_path}"', f'str(settings.paths.{self._path_to_setting(old_path)})')
            content = content.replace(f"'{old_path}'", f'str(settings.paths.{self._path_to_setting(old_path)})')
        
        return content
    
    def _path_to_setting(self, path: str) -> str:
        """Convert path to settings attribute."""
        if "datasets" in path:
            return "datasets_dir"
        elif "desktop" in path:
            return "desktop_dir"
        elif "MGCs" in path:
            return "mgc_dir"
        elif "arabidopsis" in path:
            return "arabidopsis_dir"
        else:
            return "base_data_dir"
    
    def _add_new_imports(self, content: str, target_module: str) -> str:
        """Add imports for new architecture."""
        new_imports = [
            "from ..config.settings import get_settings",
            "from ..core.base import BioinformaticsProcessor",
            "from ..core.exceptions import ProcessingError, ValidationError",
            "from ..utils.file_operations import file_manager",
            "from loguru import logger",
        ]
        
        # Add module-specific imports
        if target_module == "genomics":
            new_imports.append("from ..genomics.blast_analysis import BlastAnalyzer")
        elif target_module == "metabolic":
            new_imports.append("from ..metabolic.kegg_integration import KEGGAnalyzer")
        
        # Insert imports at the top
        lines = content.split('\n')
        import_line = -1
        
        # Find where to insert imports
        for i, line in enumerate(lines):
            if line.strip().startswith('import ') or line.strip().startswith('from '):
                import_line = i
                break
        
        if import_line >= 0:
            lines[import_line:import_line] = new_imports
        else:
            lines = new_imports + lines
        
        return '\n'.join(lines)
    
    def _transform_to_oop(self, content: str, analysis: ScriptAnalysis) -> str:
        """Transform procedural code to OOP pattern."""
        # If already has classes, return as is
        if analysis.classes:
            return content
        
        # Create class wrapper
        class_name = self._generate_class_name(analysis.script_path)
        
        oop_template = f'''
class {class_name}(BioinformaticsProcessor):
    """Migrated from legacy script: {analysis.script_path.name}."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass
'''
        
        # This is a simplified transformation
        # In practice, you'd need more sophisticated AST manipulation
        return oop_template + "\n\n" + content
    
    def _generate_class_name(self, script_path: Path) -> str:
        """Generate class name from script path."""
        name = script_path.stem
        name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
        name = ''.join(word.capitalize() for word in name.split('_'))
        if not name.endswith('Processor'):
            name += 'Processor'
        return name
    
    def _add_error_handling(self, content: str) -> str:
        """Add proper error handling."""
        # Add try-catch blocks around main operations
        # This is a simplified implementation
        return content
    
    def _add_logging(self, content: str) -> str:
        """Add logging statements."""
        # Replace print statements with logger calls
        content = re.sub(r'print\(([^)]+)\)', r'logger.info(\1)', content)
        return content
    
    def _get_target_path(self, script_path: Path, target_module: str) -> Path:
        """Get target path for migrated script."""
        base_dir = Path("src/plant_mgc_analysis")
        target_dir = base_dir / target_module
        
        # Generate new filename
        new_name = script_path.stem + "_migrated.py"
        return target_dir / new_name
    
    def migrate_batch(
        self,
        migration_plan: MigrationPlan,
        priority: str = "high",
        dry_run: bool = True
    ) -> List[Dict[str, any]]:
        """Migrate batch of scripts based on priority."""
        if priority == "high":
            scripts = migration_plan.high_priority
        elif priority == "medium":
            scripts = migration_plan.medium_priority
        elif priority == "low":
            scripts = migration_plan.low_priority
        else:
            scripts = migration_plan.analyzed_scripts
        
        results = []
        for analysis in scripts:
            result = self.migrate_script(
                analysis.script_path,
                analysis.suggested_module,
                dry_run=dry_run
            )
            results.append(result)
        
        return results
    
    def create_migration_report(self, migration_plan: MigrationPlan) -> str:
        """Create detailed migration report."""
        report = f"""
# Plant MGC Analysis Pipeline - Migration Report

## Summary
{migration_plan.get_summary()}

## High Priority Scripts ({len(migration_plan.high_priority)})
"""
        
        for analysis in migration_plan.high_priority:
            report += f"- {analysis.script_path.name} -> {analysis.suggested_module} (complexity: {analysis.complexity_score})\n"
        
        report += f"\n## Medium Priority Scripts ({len(migration_plan.medium_priority)})\n"
        for analysis in migration_plan.medium_priority:
            report += f"- {analysis.script_path.name} -> {analysis.suggested_module} (complexity: {analysis.complexity_score})\n"
        
        report += f"\n## Low Priority Scripts ({len(migration_plan.low_priority)})\n"
        for analysis in migration_plan.low_priority:
            report += f"- {analysis.script_path.name} -> {analysis.suggested_module} (complexity: {analysis.complexity_score})\n"
        
        report += f"\n## Duplicated Functions\n"
        for func, paths in migration_plan.duplications.items():
            if len(paths) > 1:
                report += f"- {func}: {len(paths)} occurrences\n"
        
        return report


def create_migration_plan(directory: Path) -> MigrationPlan:
    """Create migration plan for directory."""
    analyzer = LegacyScriptAnalyzer()
    return analyzer.analyze_directory(directory)


def migrate_legacy_scripts(
    directory: Path,
    priority: str = "high",
    dry_run: bool = True
) -> List[Dict[str, any]]:
    """Migrate legacy scripts in directory."""
    # Create migration plan
    plan = create_migration_plan(directory)
    
    # Migrate scripts
    migrator = ScriptMigrator()
    results = migrator.migrate_batch(plan, priority, dry_run)
    
    # Create report
    report = migrator.create_migration_report(plan)
    
    # Save report
    report_path = directory / "migration_report.md"
    with open(report_path, 'w') as f:
        f.write(report)
    
    return results