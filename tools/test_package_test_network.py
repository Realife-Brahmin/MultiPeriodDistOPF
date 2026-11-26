"""Packaging tests only; these do not run OPF or power-flow simulations."""
import tempfile
import unittest
from pathlib import Path
import posixpath
import re
import zipfile

from package_test_network import PROFILES, ROOT, package_network


class NetworkPackageTests(unittest.TestCase):
    def test_real_deck_and_profiles_are_byte_identical(self):
        systems = sorted(p.name for p in (ROOT / "test_networks").iterdir()
                         if p.is_dir() and (p / "Master.dss").is_file())
        self.assertIn("ieee123C_1ph", systems)
        self.assertIn("large10kC_1ph", systems)
        for system in systems:
            with self.subTest(system=system), tempfile.TemporaryDirectory() as tmp:
                output = Path(tmp) / "case.zip"
                package_network(system, output)
                with zipfile.ZipFile(output) as archive:
                    for path in (ROOT / "test_networks" / system).rglob("*"):
                        if path.is_file():
                            name = f"test_networks/{system}/{path.relative_to(ROOT / 'test_networks' / system).as_posix()}"
                            self.assertEqual(archive.read(name), path.read_bytes())
                    for profile in PROFILES:
                        self.assertEqual(archive.read(f"test_networks/{profile}"),
                                         (ROOT / "test_networks" / profile).read_bytes())
                    self.assertIn("CITATION.cff", archive.namelist())
                    self.assertIn("LICENSE", archive.namelist())
                    # Preserve redirect closure, including shared parent profiles.
                    for name in archive.namelist():
                        if not name.lower().endswith(".dss"):
                            continue
                        text = archive.read(name).decode("utf-8-sig")
                        pattern = r'(?im)^\s*(?:Redirect|Compile)\s+(?:"([^"]+)"|\[([^\]]+)\]|([^\s!]+))'
                        for match in re.finditer(pattern, text):
                            target = next(x for x in match.groups() if x is not None)
                            resolved = posixpath.normpath(posixpath.join(posixpath.dirname(name), target))
                            self.assertIn(resolved, archive.namelist(), (name, target))

    def test_deterministic_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "a.zip", Path(tmp) / "b.zip"
            package_network("ieee123C_1ph", first)
            package_network("ieee123C_1ph", second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with self.assertRaises(FileExistsError):
                package_network("ieee123C_1ph", first)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_reject_paths_missing_decks_and_input_tree_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "case.zip"
            for system in ("", ".", "..", "C:", "../ieee123C_1ph", "a/b", "a\\b", "not_a_network"):
                with self.subTest(system=system), self.assertRaises(ValueError):
                    package_network(system, output)
            with self.assertRaises(ValueError):
                package_network("ieee123C_1ph", ROOT / "test_networks" / "case.zip")
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
