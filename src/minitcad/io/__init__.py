"""输入输出：单位换算、配置读取与 long-form CSV schema。

公开 API（W2/W3 冻结口径）：

* 单位：``UnitError``、``UnitSystem``、``convert``、``cm_to_si``、``si_to_cm``、
  ``require_suffixed_field``、``CM_SUFFIXES``、``SI_SUFFIXES``
* 配置：``load_config``、``parse_config_dict``
* CSV：``CSV_COLUMNS``、``header_line``、``format_row``、``write_rows``、
  ``read_header``、``check_header``
"""

from minitcad.io.units import (
    CM,
    CM_SUFFIXES,
    SI,
    SI_SUFFIXES,
    UnitError,
    UnitSystem,
    cm_to_si,
    convert,
    require_suffixed_field,
    si_to_cm,
)
from minitcad.io.config import load_config, parse_config_dict
from minitcad.io.results_csv import (
    CSV_COLUMNS,
    RECORD_TYPES,
    REFERENCE_KINDS,
    STATUSES,
    THEORY_KINDS,
    build_rows,
    check_header,
    columns_for,
    format_row,
    header_line,
    read_header,
    write_comparison_csv,
    write_rows,
)
from minitcad.io.plot import write_comparison_png
from minitcad.io.report import validation_metrics, write_convergence_md

__all__ = [
    "CM",
    "SI",
    "CM_SUFFIXES",
    "SI_SUFFIXES",
    "UnitError",
    "UnitSystem",
    "convert",
    "cm_to_si",
    "si_to_cm",
    "require_suffixed_field",
    "load_config",
    "parse_config_dict",
    "CSV_COLUMNS",
    "RECORD_TYPES",
    "THEORY_KINDS",
    "REFERENCE_KINDS",
    "STATUSES",
    "header_line",
    "format_row",
    "build_rows",
    "columns_for",
    "write_rows",
    "write_comparison_csv",
    "write_comparison_png",
    "write_convergence_md",
    "validation_metrics",
    "read_header",
    "check_header",
]
