import yaml
import os


def load_paths():
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "config",
        "paths.yaml"
    )

    with open(config_path, "r") as f:
        return yaml.safe_load(f)
