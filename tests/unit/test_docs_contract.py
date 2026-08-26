"""Contracts that keep user-facing documentation aligned with the CLI."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shlex
import unittest


from bsreadsim.cli import build_parser, build_run_document


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPOSITORY_ROOT / "docs"


def documented_bash_commands() -> list[tuple[Path, str]]:
    """Return complete continued ``bsreadsim`` commands from Markdown blocks."""
    markdown_files = [REPOSITORY_ROOT / "README.md"]
    markdown_files.extend(sorted(DOCS_ROOT.rglob("*.md")))
    command_pattern = re.compile(r"(?m)^bsreadsim(?:[^\n]*\\\n)*[^\n]*")
    commands: list[tuple[Path, str]] = []
    for path in markdown_files:
        text = path.read_text(encoding="utf-8")
        for block in re.findall(r"```bash\n(.*?)```", text, re.DOTALL):
            commands.extend(
                (path, command.replace("\\\n", " "))
                for command in command_pattern.findall(block)
            )
    return commands


def public_long_options(parser: argparse.ArgumentParser) -> set[str]:
    """Collect public long options from a parser and all of its subparsers."""
    options: set[str] = set()
    for action in parser._actions:
        options.update(
            option for option in action.option_strings if option.startswith("--")
        )
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                options.update(public_long_options(subparser))
    return options


def public_short_options(parser: argparse.ArgumentParser) -> set[str]:
    """Collect public one-letter aliases from a parser and its subparsers."""
    options: set[str] = set()
    for action in parser._actions:
        options.update(
            option
            for option in action.option_strings
            if re.fullmatch(r"-[A-Za-z]", option)
        )
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                options.update(public_short_options(subparser))
    return options


class DocumentationContractTests(unittest.TestCase):
    def test_primary_navigation_keeps_the_six_entry_points(self) -> None:
        navigation = (REPOSITORY_ROOT / "mkdocs.yml").read_text(encoding="utf-8")
        entries = (
            "  - Overview: index.md",
            "  - Installation: getting-started/installation.md",
            "  - Quick start: getting-started/quickstart.md",
            "  - Simulation:",
            "      - Workflow: simulation/workflow.md",
            "      - Tutorials: simulation/tutorials.md",
            "      - Customize: simulation/customize.md",
            "      - Other assays: simulation/other-assays.md",
            "  - Outputs:",
            "      - Overview: outputs/index.md",
            "      - Inspect outputs: outputs/inspect.md",
            "  - Reference:",
        )
        positions = [navigation.find(entry) for entry in entries]

        self.assertNotIn(-1, positions)
        self.assertEqual(positions, sorted(positions))

    def test_simulation_section_separates_workflow_configuration_and_tutorials(self) -> None:
        workflow = (DOCS_ROOT / "simulation" / "workflow.md").read_text(
            encoding="utf-8"
        )
        configuration = (DOCS_ROOT / "simulation" / "customize.md").read_text(
            encoding="utf-8"
        )
        tutorials = (DOCS_ROOT / "simulation" / "tutorials.md").read_text(
            encoding="utf-8"
        )
        other_assays = (DOCS_ROOT / "simulation" / "other-assays.md").read_text(
            encoding="utf-8"
        )

        workflow_headings = (
            "# Workflow",
            "## From reference genome to sequencing reads",
            "## Five stages of read simulation",
            "## Six supported sequencing assays",
        )
        configuration_headings = (
            "# Customize",
            "## Genome { #genetic-variation }",
            "## Methylome { #methylation }",
            "### Probability sources",
            "### Methylation pattern",
            "### Allele-specific methylation",
            "## Fragmentation { #supported-technologies }",
            "### WGBS { #wgbs }",
            "### RRBS { #rrbs }",
            "### TBS { #tbs }",
            "## Sequencing { #library-and-sequencing }",
            "## Output { #reproducibility }",
        )
        tutorial_headings = (
            "# Tutorials",
            "## Synthetic WGBS { #wgbs }",
            "## Variant sets and methylation profiles",
            "## Allele-specific methylation",
            "## RRBS and TBS",
            "### RRBS { #rrbs }",
            "### TBS { #tbs }",
            "## Ground-truth benchmarking",
        )
        other_assay_headings = (
            "# Other assays",
            "## How standard sequencing is simulated",
            "## Choose an assay",
            "### WGS { #wgs }",
            "### WES { #wes }",
            "### TS { #ts }",
            "## Shared simulation controls",
            "## Outputs and ground truth",
        )

        for page, headings in (
            (workflow, workflow_headings),
            (configuration, configuration_headings),
            (tutorials, tutorial_headings),
            (other_assays, other_assay_headings),
        ):
            positions = [page.find(heading) for heading in headings]
            self.assertNotIn(-1, positions)
            self.assertEqual(positions, sorted(positions))

        self.assertIn('class="workflow-figure"', workflow)
        self.assertIn("../../img/BSReadSim_workflow.png", workflow)
        self.assertIn("../../img/BS_seqtech.png", workflow)
        self.assertNotIn("## 1. Introduce genetic variation", workflow)
        self.assertNotIn("### WGS { #wgs }", configuration)
        self.assertNotIn("### WES { #wes }", configuration)
        self.assertNotIn("### TS { #ts }", configuration)
        self.assertEqual(
            sorted(path.name for path in (DOCS_ROOT / "simulation").glob("*.md")),
            ["customize.md", "other-assays.md", "tutorials.md", "workflow.md"],
        )

    def test_quickstart_routes_real_runs_to_simulation_sections(self) -> None:
        quickstart = (DOCS_ROOT / "getting-started" / "quickstart.md").read_text(
            encoding="utf-8"
        )

        self.assertIn("../simulation/customize.md#wgbs", quickstart)
        self.assertIn("../simulation/customize.md#rrbs", quickstart)
        self.assertIn("../simulation/customize.md#tbs", quickstart)
        self.assertIn("../simulation/other-assays.md", quickstart)
        self.assertNotIn("../simulation/other-assays.md#wgs", quickstart)
        self.assertNotIn("../simulation/other-assays.md#wes", quickstart)
        self.assertNotIn("../simulation/other-assays.md#ts", quickstart)
        self.assertIn("../simulation/customize.md", quickstart)
        self.assertIn("../simulation/tutorials.md", quickstart)

    def test_command_line_overview_mentions_every_public_option(self) -> None:
        overview = (DOCS_ROOT / "reference" / "cli.md").read_text(encoding="utf-8")
        missing = sorted(
            option for option in public_long_options(build_parser()) if option not in overview
        )

        self.assertEqual(missing, [])

        missing_short = sorted(
            option
            for option in public_short_options(build_parser())
            if "`{}`".format(option) not in overview
        )
        self.assertEqual(missing_short, [])

    def test_quickstart_command_projects_to_the_default_wgbs_model(self) -> None:
        quickstart = (DOCS_ROOT / "getting-started" / "quickstart.md").read_text(
            encoding="utf-8"
        )
        command_match = re.search(
            r"```bash\n(bsreadsim run wgbs [^\n]+)\n```", quickstart
        )
        self.assertIsNotNone(command_match)
        command = command_match.group(1)
        argv = shlex.split(command)

        arguments = build_parser().parse_args(argv[1:])
        document = build_run_document(arguments, REPOSITORY_ROOT)

        self.assertEqual(document["technology"], "WGBS")
        self.assertEqual(document["mutation"]["rate"], 0.001)
        self.assertEqual(document["fragments"]["count"], 1000)

    def test_documented_run_commands_choose_a_variant_policy(self) -> None:
        missing = []
        for path, command in documented_bash_commands():
            if "bsreadsim run " not in command or "--help" in command:
                continue
            if path == DOCS_ROOT / "getting-started" / "quickstart.md":
                continue
            if "--mutation-rate" not in command and "--vcf" not in command:
                missing.append("{}: {}".format(path, command))

        self.assertEqual(missing, [])

    def test_documented_cli_commands_parse_with_the_current_interface(self) -> None:
        parser = build_parser()
        failures = []
        for path, command in documented_bash_commands():
            if "--help" in command or "--version" in command:
                continue
            try:
                parser.parse_args(shlex.split(command)[1:])
            except (SystemExit, ValueError) as error:
                failures.append("{}: {} ({})".format(path, command, error))

        self.assertEqual(failures, [])

    def test_removed_insert_size_option_is_absent_from_user_docs(self) -> None:
        markdown_files = [REPOSITORY_ROOT / "README.md"]
        markdown_files.extend(sorted(DOCS_ROOT.rglob("*.md")))
        stale = [
            str(path.relative_to(REPOSITORY_ROOT))
            for path in markdown_files
            if "--insert-size" in path.read_text(encoding="utf-8")
        ]

        self.assertEqual(stale, [])

    def test_user_docs_use_consistent_prepared_truth_terms(self) -> None:
        markdown_files = [REPOSITORY_ROOT / "README.md"]
        markdown_files.extend(sorted(DOCS_ROOT.rglob("*.md")))
        stale_patterns = (
            re.compile(r"\bvariant catalog\b", re.IGNORECASE),
            re.compile(r"\bmethylation catalog\b", re.IGNORECASE),
            re.compile(r"\bcandidate catalog\b", re.IGNORECASE),
            re.compile(r"\bbsreadsim\s+catalog\b", re.IGNORECASE),
        )
        stale = []
        for path in markdown_files:
            text = path.read_text(encoding="utf-8")
            for pattern in stale_patterns:
                if pattern.search(text):
                    stale.append(
                        "{}: {}".format(
                            path.relative_to(REPOSITORY_ROOT), pattern.pattern
                        )
                    )

        self.assertEqual(stale, [])

        formats = (DOCS_ROOT / "reference" / "formats.md").read_text(
            encoding="utf-8"
        )
        outputs = (DOCS_ROOT / "outputs" / "index.md").read_text(
            encoding="utf-8"
        )
        for term in ("prepared variant set", "prepared methylation profile"):
            self.assertIn(term, formats.lower())
        for term in ("saved variant set", "saved methylation profile"):
            self.assertIn(term, outputs.lower())
        self.assertIn("simulation truth artifacts", outputs.lower())

    def test_depth_documentation_uses_the_frozen_ceil_contract(self) -> None:
        cli = (DOCS_ROOT / "reference" / "cli.md").read_text(encoding="utf-8")
        configuration = (DOCS_ROOT / "simulation" / "customize.md").read_text(
            encoding="utf-8"
        )

        self.assertIn("ceil(effective_reference_bases", cli)
        self.assertIn("ceil(effective_reference_bases", configuration)

    def test_quickstart_reference_supports_the_default_insert_range(self) -> None:
        quickstart = (DOCS_ROOT / "getting-started" / "quickstart.md").read_text(
            encoding="utf-8"
        )

        self.assertIn(
            "bsreadsim export test-fasta -o test.fa",
            quickstart,
        )
        self.assertNotIn("--insert-size", quickstart)
        command_match = re.search(
            r"```bash\n(bsreadsim run wgbs [^\n]+)\n```", quickstart
        )
        self.assertIsNotNone(command_match)
        command = command_match.group(1)
        arguments = build_parser().parse_args(shlex.split(command)[1:])
        document = build_run_document(arguments, REPOSITORY_ROOT)
        fasta = REPOSITORY_ROOT / "data" / "examples" / "test.fa"
        reference_length = sum(
            len(line.strip())
            for line in fasta.read_text(encoding="ascii").splitlines()
            if not line.startswith(">")
        )
        self.assertGreater(reference_length, document["fragments"]["insert_max"])

    def test_outputs_show_the_shape_of_every_published_artifact(self) -> None:
        outputs = (DOCS_ROOT / "outputs" / "index.md").read_text(
            encoding="utf-8"
        )
        examples = (
            "sim.R1.fastq.gz",
            "sim.R2.fastq.gz",
            "sim.bam",
            "sim.manifest.json",
            "truth/sim.methdb",
            "@chr1:101-108:0/1",
            "@chr1:101-108:0/2",
            '"status": "complete"',
            "samtools view",
            "zt:Z:AAAAA",
            "00000000: 6d65 7468 6462",
            "truth/sim.variants.vcf.gz",
        )

        self.assertEqual(
            [example for example in examples if example not in outputs],
            [],
        )


if __name__ == "__main__":
    unittest.main()
