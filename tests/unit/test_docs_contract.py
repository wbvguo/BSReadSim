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
    def test_overview_keeps_complete_html_sections(self) -> None:
        overview = (DOCS_ROOT / "index.md").read_text(encoding="utf-8")

        for fragment in (
            '<div class="docs-hero" markdown>',
            '<img class="docs-hero__logo"',
            "# BSReadSim",
            '<div class="docs-actions" markdown>',
            "## Supported assays",
            '<div class="technology-grid" markdown>',
            "## Citation",
            "@article{guo2024bsreadsim,",
        ):
            self.assertIn(fragment, overview)
        self.assertEqual(overview.count("<div"), overview.count("</div>"))
        self.assertEqual(overview.count('<div class="technology-card'), 4)

    def test_customize_uses_the_cli_option_table_schema(self) -> None:
        lines = (DOCS_ROOT / "simulation" / "customize.md").read_text(
            encoding="utf-8"
        ).splitlines()
        header = "| Option | Value | Default | Description |"
        separator = "| --- | --- | --- | --- |"
        table_starts = [index for index, line in enumerate(lines) if line == header]

        self.assertGreater(len(table_starts), 0)
        self.assertNotIn("| Option | Default | Effect |", lines)
        for start in table_starts:
            self.assertEqual(lines[start + 1], separator)
            row = start + 2
            while row < len(lines) and lines[row].startswith("|"):
                self.assertEqual(
                    len(lines[row].strip("|").split("|")),
                    4,
                    msg="Customize option row {} does not have four columns".format(
                        row + 1
                    ),
                )
                row += 1

    def test_flags_follow_valued_options_in_option_tables(self) -> None:
        for path in (
            DOCS_ROOT / "simulation" / "customize.md",
            DOCS_ROOT / "reference" / "cli.md",
        ):
            lines = path.read_text(encoding="utf-8").splitlines()
            for start, line in enumerate(lines):
                if line != "| Option | Value | Default | Description |":
                    continue
                seen_flag = False
                row = start + 2
                while row < len(lines) and lines[row].startswith("|"):
                    columns = lines[row].strip("|").split("|")
                    value = columns[1].strip()
                    if value == "Flag":
                        seen_flag = True
                    else:
                        self.assertFalse(
                            seen_flag,
                            msg="{}:{} places a valued option after a flag".format(
                                path.relative_to(REPOSITORY_ROOT), row + 1
                            ),
                        )
                    row += 1

    def test_primary_navigation_keeps_entry_points_in_order(self) -> None:
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
            "  - Outputs: outputs/index.md",
            "  - Reference:",
            "      - Troubleshoot: help/troubleshoot.md",
        )
        positions = [navigation.find(entry) for entry in entries]

        self.assertNotIn(-1, positions)
        self.assertEqual(positions, sorted(positions))
        self.assertIn("    - navigation.expand", navigation)
        self.assertNotIn("    - navigation.sections", navigation)
        self.assertNotIn("getting-started/platforms.md", navigation)

    def test_installation_owns_the_platform_support_contract(self) -> None:
        installation = (DOCS_ROOT / "getting-started" / "installation.md").read_text(
            encoding="utf-8"
        )

        self.assertIn("## Supported platforms", installation)
        for environment in (
            "Linux x86-64",
            "Windows through WSL2",
            "Native Windows",
            "macOS",
            "Linux ARM64",
        ):
            self.assertIn(environment, installation)
        self.assertFalse((DOCS_ROOT / "getting-started" / "platforms.md").exists())

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
            "### Generate de novo variants { #generated-variants }",
            "### Load variants from a VCF { #vcf-genome }",
            "## Methylation { #methylation }",
            "### Generate methylation probabilities { #generated-methylation }",
            "### Load methylation profiles { #predefined-methylation }",
            "### Add allele-specific methylation { #allele-specific-methylation }",
            "### Realize methylation states { #methylation-states }",
            "## Fragment sampling { #supported-technologies }",
            "### Generate fragments",
            "### Set fragment length { #fragment-geometry }",
            "### Sample fragments",
            "## Bisulfite conversion { #bisulfite-conversion }",
            "## Read generation { #library-and-sequencing }",
            "### Set number of reads { #dataset-size }",
            "### Set read layout { #read-layout }",
            "### Configure base quality and sequencing errors { #quality-and-error }",
            "## Output and reproducibility { #reproducibility }",
            "### Configure execution and output { #output-format }",
            "### Save simulation truth { #truth-artifacts }",
            "### Control random variation { #random-seeds }",
        )
        tutorial_headings = (
            "# Tutorials",
            "## WGBS from a reference genome { #wgbs }",
            "## Variant sets and methylation profiles",
            "## Allele-specific methylation",
            "## Enrichment-based assays",
            "### RRBS { #rrbs }",
            "### TBS { #tbs }",
            "## Ground truth in simulation",
            "### Save",
            "### Reuse",
        )
        other_assay_headings = (
            "# Other assays",
            "## Non-bisulfite workflow { #how-standard-sequencing-is-simulated }",
            "## Choose an assay",
            "### WGS { #wgs }",
            "### WES { #wes }",
            "### TS { #ts }",
            "## Shared simulation controls",
            "## Outputs and simulation truth",
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

    def test_cli_reference_separates_rules_from_core_options(self) -> None:
        overview = (DOCS_ROOT / "reference" / "cli.md").read_text(encoding="utf-8")

        self.assertIn(
            "## Option combination rules { #option-combination-rules }",
            overview,
        )
        self.assertIn(
            "## Core integration options { #core-integration-options }",
            overview,
        )
        self.assertNotIn("## Compatibility", overview)
        for retired_alias in (
            "--bed-methyl",
            "--pool-methylation-values",
            "--cgmap-pool",
        ):
            self.assertNotIn(retired_alias, overview)

    def test_quickstart_command_projects_to_the_default_wgbs_model(self) -> None:
        quickstart = (DOCS_ROOT / "getting-started" / "quickstart.md").read_text(
            encoding="utf-8"
        )
        command_match = re.search(
            r"```bash\n(bsreadsim run wgbs(?:[^\n]*\\\n)*[^\n]*)\n```",
            quickstart,
        )
        self.assertIsNotNone(command_match)
        command = command_match.group(1).replace("\\\n", " ")
        argv = shlex.split(command)

        arguments = build_parser().parse_args(argv[1:])
        document = build_run_document(arguments, REPOSITORY_ROOT)

        self.assertEqual(document["technology"], "WGBS")
        self.assertEqual(document["mutation"]["rate"], 0.001)
        self.assertEqual(document["reads"]["count"], 1000)
        self.assertEqual(document["fragments"]["read_length_1"], 100)
        self.assertEqual(document["fragments"]["insert_mean"], 400)
        self.assertEqual(document["fragments"]["insert_sd"], 25.0)
        self.assertEqual(document["sequencing"]["conversion_rate"], 0.998)
        self.assertEqual(document["execution"]["threads"], 4)
        for option in (
            "`run wgbs`",
            "`-r`",
            "`-o`",
            "`-n`",
            "`-t`",
            "`--seed`",
        ):
            self.assertIn(option, quickstart)

    def test_inline_long_option_references_are_public(self) -> None:
        parser_options = public_long_options(build_parser())
        markdown_files = [REPOSITORY_ROOT / "README.md"]
        markdown_files.extend(sorted(DOCS_ROOT.rglob("*.md")))
        stale: list[str] = []
        pattern = re.compile(r"`(--[a-z][a-z0-9-]*)`")
        for path in markdown_files:
            for option in pattern.findall(path.read_text(encoding="utf-8")):
                if option not in parser_options:
                    stale.append(
                        "{}: {}".format(path.relative_to(REPOSITORY_ROOT), option)
                    )

        self.assertEqual(stale, [])

    def test_tbs_uses_the_shared_variable_insert_default(self) -> None:
        arguments = build_parser().parse_args(
            shlex.split(
                "bsreadsim run tbs -r reference.fa -o runs/tbs -n 2 "
                "--targets targets.bed --mutation-rate 0"
            )[1:]
        )
        document = build_run_document(arguments, REPOSITORY_ROOT)
        troubleshoot = (DOCS_ROOT / "help" / "troubleshoot.md").read_text(
            encoding="utf-8"
        )

        self.assertEqual(document["fragments"]["insert_mean"], 400)
        self.assertEqual(document["fragments"]["insert_sd"], 25.0)
        self.assertIn("For WGBS, WGS, TBS, WES, and TS", troubleshoot)
        self.assertIn("With a positive\n  SD", troubleshoot)

    def test_troubleshoot_is_symptom_oriented_and_uses_current_terms(self) -> None:
        path = DOCS_ROOT / "help" / "troubleshoot.md"
        troubleshoot = path.read_text(encoding="utf-8")

        self.assertTrue(troubleshoot.startswith("# Troubleshoot\n"))
        self.assertFalse((DOCS_ROOT / "help" / "troubleshooting.md").exists())
        headings = (
            "## Commands and setup",
            "## Input files",
            "## Read generation and sampling",
            "## Outputs",
            "## Performance",
        )
        positions = [troubleshoot.find(heading) for heading in headings]
        self.assertNotIn(-1, positions)
        self.assertEqual(positions, sorted(positions))
        for current_rule in (
            "### `bsreadsim` is not found",
            "bsreadsim validate \\",
            "only a subset of FASTA contigs",
            "Use `--strict`",
            "`--sampling gc` and `--gc-profile PATH`",
            (
                "VCF, ASM, nonzero de novo mutation, or a MethDB containing "
                "embedded variants"
            ),
            "drops unreachable\npositive mass and renormalizes",
            "do not need to match the later run",
            "the parent directory must already exist",
            "has enough free\nspace",
            "Compare `-t 1`, `-t 2`, and `-t 4`",
        ):
            self.assertIn(current_rule, troubleshoot)
        self.assertNotIn("the fixed fragment and read lengths", troubleshoot)

    def test_documented_run_commands_choose_a_variant_policy(self) -> None:
        missing = []
        for path, command in documented_bash_commands():
            if "bsreadsim run " not in command or "--help" in command:
                continue
            if path == DOCS_ROOT / "getting-started" / "quickstart.md":
                continue
            if not any(
                policy in command
                for policy in (
                    "--mutation-rate",
                    "--vcf",
                    "--methdb",
                    "--asm",
                    "--asm-bed",
                )
            ):
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

    def test_removed_fragments_option_is_absent_from_user_docs(self) -> None:
        markdown_files = [REPOSITORY_ROOT / "README.md"]
        markdown_files.extend(sorted(DOCS_ROOT.rglob("*.md")))
        stale = [
            str(path.relative_to(REPOSITORY_ROOT))
            for path in markdown_files
            if "--fragments" in path.read_text(encoding="utf-8")
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
        for term in ("exported vcf", "exported methdb"):
            self.assertIn(term, outputs.lower())

    def test_file_format_contracts_live_on_one_central_page(self) -> None:
        formats = (DOCS_ROOT / "reference" / "formats.md").read_text(
            encoding="utf-8"
        )
        headings = (
            "# File formats",
            "## Genome and variants { #genome-and-variant-inputs }",
            "## Methylation profiles { #methylation-profile-inputs }",
            "## Allele-specific methylation { #allele-specific-methylation-inputs }",
            "## Fragment generation and sampling { #fragment-sampling-inputs }",
            "## Sequencing models { #sequencing-model-inputs }",
        )

        positions = [formats.find(heading) for heading in headings]
        self.assertNotIn(-1, positions)
        self.assertEqual(positions, sorted(positions))
        self.assertEqual(
            sorted((DOCS_ROOT / "reference" / "formats").glob("*.md")),
            [],
        )
        self.assertNotIn("## Generated files", formats)

    def test_file_format_contracts_cover_parser_boundaries(self) -> None:
        formats = (DOCS_ROOT / "reference" / "formats.md").read_text(
            encoding="utf-8"
        )

        for text in (
            "| RRBS candidate BED | Plain text only |",
            "| Quality and error models | Uncompressed JSON only |",
            "For input files, the",
            "Plain and gzip-compressed FASTA are",
            "contig length representable as uint32",
            "A header-only VCF is valid\nand represents an empty variant set.",
            "they are skipped before exact\nREF matching",
            "repeated IDs are\ndisambiguated",
            "Every row must use the same supported width: 6, 8, 9, or 10 fields.",
            "#chrom\tchromStart\tchromEnd\tprobability",
            "#CHR\tNUC\tPOS\tCONTEXT\tDINUC\tMETH\tMC\tNC",
            "Only `percentModified` determines",
            "`METH` is the authoritative",
            "not the `--methbed` input format",
            "Exactly one of `Allele1` and `Allele2`",
            "none accepts\n`.` or `na`",
            "For fixed-insert runs, positive probability",
            "variable-insert runs instead drop unreachable",
            "Rows may appear in any order;",
            "Every count is uint32",
            "does not accept MethDB",
        ):
            self.assertIn(text, formats)

        self.assertEqual(
            formats.count("| Column | Name | Requirement |"),
            8,
        )

        self.assertNotIn(
            "A requested bin with no eligible fragments is rejected.",
            formats,
        )

    def test_external_input_formats_link_to_upstream_definitions(self) -> None:
        formats = (DOCS_ROOT / "reference" / "formats.md").read_text(
            encoding="utf-8"
        )

        for upstream_url in (
            "https://www.ncbi.nlm.nih.gov/genbank/fastaformat/",
            "https://samtools.github.io/hts-specs/VCFv4.2.pdf",
            "https://samtools.github.io/hts-specs/VCFv4.3.pdf",
            "https://genome.ucsc.edu/goldenPath/help/bedMethyl.html",
            "https://cgmaptools.github.io/cgmaptools_documentation/file-formats.html#cgmap-format",
            "https://cgmaptools.github.io/cgmaptools_documentation/methylation-analysis.html#asm",
            "https://genome.ucsc.edu/FAQ/FAQformat.html#format1",
        ):
            self.assertIn(upstream_url, formats)

        self.assertIn("### ASM { #asm }", formats)
        self.assertNotIn("### CGmapTools ASM", formats)

    def test_depth_documentation_uses_the_reads_first_contract(self) -> None:
        cli = (DOCS_ROOT / "reference" / "cli.md").read_text(encoding="utf-8")
        configuration = (DOCS_ROOT / "simulation" / "customize.md").read_text(
            encoding="utf-8"
        )

        formulas = (
            "raw_reads = effective_reference_bases * D / read_length",
            "resolved_reads = emitted_mates * ceil(raw_reads / emitted_mates)",
            "resolved_fragments = resolved_reads / emitted_mates",
        )
        for formula in formulas:
            self.assertIn(
                formula,
                cli,
            )
            self.assertNotIn(formula, configuration)
        self.assertIn(
            "../reference/cli.md#required-inputs-and-dataset-size",
            configuration,
        )

    def test_quickstart_reference_supports_the_default_insert_range(self) -> None:
        quickstart = (DOCS_ROOT / "getting-started" / "quickstart.md").read_text(
            encoding="utf-8"
        )

        self.assertNotIn("bsreadsim export test-fasta", quickstart)
        self.assertIn(
            "https://github.com/wbvguo/BSReadSim/blob/main/data/examples/test.fa",
            quickstart,
        )
        self.assertNotIn("--insert-size", quickstart)
        command_match = re.search(
            r"```bash\n(bsreadsim run wgbs(?:[^\n]*\\\n)*[^\n]*)\n```",
            quickstart,
        )
        self.assertIsNotNone(command_match)
        command = command_match.group(1).replace("\\\n", " ")
        arguments = build_parser().parse_args(shlex.split(command)[1:])
        document = build_run_document(arguments, REPOSITORY_ROOT)
        fasta = REPOSITORY_ROOT / "data" / "examples" / "test.fa"
        reference_length = sum(
            len(line.strip())
            for line in fasta.read_text(encoding="ascii").splitlines()
            if not line.startswith(">")
        )
        self.assertGreater(reference_length, document["fragments"]["insert_max"])

    def test_output_details_live_on_one_outputs_page(self) -> None:
        outputs = (DOCS_ROOT / "outputs" / "index.md").read_text(
            encoding="utf-8"
        )
        formats = (DOCS_ROOT / "reference" / "formats.md").read_text(
            encoding="utf-8"
        )
        examples = (
            "sim.R1.fastq.gz",
            "sim.R2.fastq.gz",
            "sim.bam",
            "sim.manifest.json",
            "sim.methdb",
            '"status": "complete"',
            "sim.variants.vcf.gz",
        )

        self.assertEqual(
            [example for example in examples if example not in outputs],
            [],
        )
        self.assertEqual(
            sorted(path.name for path in (DOCS_ROOT / "outputs").glob("*.md")),
            ["index.md"],
        )
        self.assertNotIn("## Choose other outputs", outputs)
        for text in ("WGBS, RRBS, and TBS", "python -m json.tool", "samtools view"):
            self.assertIn(text, outputs)
        self.assertIn("@chr1:101-108:2a/1", outputs)
        for text in ("## Run manifest", "## Truth artifacts"):
            self.assertIn(text, outputs)
        for tag in ("XG:Z", "XR:Z", "YS:Z", "zt:Z"):
            self.assertIn(tag, outputs)
        self.assertFalse((DOCS_ROOT / "reference" / "output-formats.md").exists())
        self.assertIn("`methdb` magic", formats)
        self.assertIn("version byte `2`", formats)
        self.assertNotIn("00000000: 6d65 7468 6462 02", outputs)


if __name__ == "__main__":
    unittest.main()
