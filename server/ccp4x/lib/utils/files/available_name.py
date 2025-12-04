import pathlib


def available_file_name_based_on(file_path: pathlib.Path):
    destination_dir = file_path.parent
    base_stem = file_path.stem
    base_suffix = file_path.suffix
    print(file_path, destination_dir, base_stem, base_suffix)
    dest = (destination_dir / base_stem).with_suffix(base_suffix).absolute()
    attempt_number = 0
    while dest.exists():
        attempt_number += 1
        dest = (
            (destination_dir / f"{base_stem}_{attempt_number}")
            .with_suffix(base_suffix)
            .absolute()
        )
    # Check that the filename does lie in the appropriate direction
    return dest
