"""Invalid user inputs should produce useful errors without a native solve."""


import pytest

from polyfempy import solve


def test_config_path_errors_name_the_file(tmp_path):
    missing = tmp_path / "missing.json"
    with pytest.raises(FileNotFoundError) as error:
        solve(cfg=missing)
    assert str(missing) in str(error.value)
    bad = tmp_path / "bad.json"
    bad.write_text('{"geometry": }', encoding="utf-8")
    with pytest.raises(ValueError, match="line 1 column") as error:
        solve(cfg=bad)
