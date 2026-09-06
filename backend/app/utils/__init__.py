# # app/utils/__init__.py
# """
# Exports all utility functions for app.utils
# so callers can do:
#     from app.utils import compute_sha256, update_results_json, ...
# """

# # IO helpers
# from .io_helpers import (
#     write_json,
#     read_json,
#     write_yaml,
#     read_yaml,
#     compute_sha256,
# )


# # Run directory helpers
# from .run_dir_helper import (
#     get_base_workdir,
#     ensure_run_dir,
#     ensure_step_dir,
#     update_results_json,  # <- make sure this is included
# )

# # Run ID helpers
# from .run_id import new_run_id, sanitize_run_id

# # Subprocess wrappers
# from .subprocess_wrappers import run_cmd

# __all__ = [
#     # io_helpers
#     "write_json",
#     "write_yaml",
#     "read_json",
#     "read_yaml",
#     "compute_sha256",
#     # run_dir_helper
#     "get_base_workdir",
#     "ensure_run_dir",
#     "ensure_step_dir",
#     "update_results_json",
#     # run_id
#     "new_run_id",
#     "sanitize_run_id",
#     # subprocess
#     "run_cmd",
# ]
