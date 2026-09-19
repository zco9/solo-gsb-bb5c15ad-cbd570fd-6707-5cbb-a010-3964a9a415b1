"""经验式 / 分子式反推。

由元素质量百分比出发：先按内置原子量表换算摩尔比（除以最小摩尔值归一），
再在 1～12 倍尺度内把摩尔比四舍五入为正整数下标，仅保留互质（gcd == 1）
候选作为经验式；用候选式反算各元素理论质量百分比，取最大绝对偏差不超过
请求 tolerance 的候选，偏差相同时依次按总原子数较少、Hill 顺序字符串较小
择优。提供实测摩尔质量时，取与其最接近的经验式摩尔质量正整数倍确定分子式，
相对偏差不得超过 tolerance%。全程 Decimal 高精度计算，输出统一四舍五入
（ROUND_HALF_UP）为定长小数字符串。
"""
from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP, localcontext
from math import gcd
from typing import Any

from .stoichiometry import ATOMIC_WEIGHTS, _round

# 摩尔比放大尺度范围（含端点）
SCALE_MIN = 1
SCALE_MAX = 12

# 输出小数字段统一保留 6 位小数
OUTPUT_PRECISION = 6

_SUPPORTED_ELEMENTS = "、".join(ATOMIC_WEIGHTS)
_ONE = Decimal(1)
_HUNDRED = Decimal(100)


class FormulaInferInputError(ValueError):
    """请求输入不合法（对应 HTTP 400）。"""


class NoMatchingFormulaError(Exception):
    """容差内无合格候选，或摩尔质量倍数相对偏差超限（对应 HTTP 422）。"""

    def __init__(self, message: str, **detail: Any) -> None:
        super().__init__(message)
        self.message = message
        self.detail = detail


def _hill_order(elements: list[str]) -> list[str]:
    """Hill 系统排序：含碳时 C、H 在前，其余按符号字母序；不含碳时全部字母序。"""
    if "C" in elements:
        head = ["C"] + (["H"] if "H" in elements else [])
        return head + sorted(e for e in elements if e not in ("C", "H"))
    return sorted(elements)


def _format_formula(subscripts: dict[str, int], order: list[str]) -> str:
    """按给定元素顺序拼接化学式，下标 1 省略。"""
    return "".join(
        e + ("" if subscripts[e] == 1 else str(subscripts[e])) for e in order
    )


def _gcd_all(values: list[int]) -> int:
    g = 0
    for v in values:
        g = gcd(g, v)
    return abs(g)


def _build_candidate(
    scale: int,
    ratios: dict[str, Decimal],
    given: dict[str, Decimal],
) -> dict[str, Any] | None:
    """构造某一尺度下的候选；四舍五入后下标不互质则返回 None。"""
    subscripts = {
        el: int((ratio * scale).quantize(_ONE, rounding=ROUND_HALF_UP))
        for el, ratio in ratios.items()
    }
    # 归一后最小摩尔比恰为 1，故所有下标均为正；互质才是经验式
    if _gcd_all(list(subscripts.values())) != 1:
        return None

    masses = {
        el: Decimal(n) * ATOMIC_WEIGHTS[el] for el, n in subscripts.items()
    }
    total_mass = sum(masses.values(), Decimal(0))
    theoretical = {el: masses[el] / total_mass * _HUNDRED for el in subscripts}
    # 带符号偏差：理论值 - 实测值（百分点）
    deviations = {el: theoretical[el] - given[el] for el in subscripts}
    max_abs = max(abs(v) for v in deviations.values())

    order = _hill_order(list(subscripts))
    return {
        "scale": scale,
        "subscripts": subscripts,
        "theoretical": theoretical,
        "deviations": deviations,
        "max_abs_deviation": max_abs,
        "total_atoms": sum(subscripts.values()),
        "hill_formula": _format_formula(subscripts, order),
        "hill_order": order,
        "molar_mass": total_mass,
    }


def _select_candidate(
    candidates: list[dict[str, Any]], tolerance: Decimal
) -> dict[str, Any]:
    """在合格候选中按 (最大绝对偏差, 总原子数, Hill 字符串) 升序择优。"""
    eligible = [c for c in candidates if c["max_abs_deviation"] <= tolerance]
    if not eligible:
        closest = min(candidates, key=lambda c: c["max_abs_deviation"])
        raise NoMatchingFormulaError(
            "1～12 倍尺度内不存在最大绝对偏差不超过 tolerance 的经验式候选",
            reason="no_empirical_candidate",
            tolerance=_round(tolerance, OUTPUT_PRECISION),
            best_formula=closest["hill_formula"],
            best_max_abs_deviation=_round(
                closest["max_abs_deviation"], OUTPUT_PRECISION
            ),
        )
    return min(
        eligible,
        key=lambda c: (
            c["max_abs_deviation"],
            c["total_atoms"],
            c["hill_formula"],
        ),
    )


def infer_formula(
    percentages: dict[str, Decimal],
    tolerance: Decimal,
    measured_molar_mass: Decimal | None,
) -> dict[str, Any]:
    """经验式反推主流程，返回可直接用于响应模型的字典。

    FormulaInferInputError / NoMatchingFormulaError 由路由层映射为 400 / 422。
    """
    # ---- 输入校验 ----
    if not percentages:
        raise FormulaInferInputError("elements 至少需要一种元素的质量百分比")
    unknown = [el for el in percentages if el not in ATOMIC_WEIGHTS]
    if unknown:
        raise FormulaInferInputError(
            f"元素 {'、'.join(unknown)} 不在内置原子量表"
            f"（支持：{_SUPPORTED_ELEMENTS}）"
        )
    for el, pct in percentages.items():
        if not pct.is_finite() or pct <= 0:
            raise FormulaInferInputError(
                f"元素 {el} 的质量百分比必须为正数，收到 {pct}"
            )
    if tolerance is None or not tolerance.is_finite() or tolerance <= 0:
        raise FormulaInferInputError("tolerance 必须为正数")
    if measured_molar_mass is not None and (
        not measured_molar_mass.is_finite() or measured_molar_mass <= 0
    ):
        raise FormulaInferInputError(
            f"molar_mass 必须为正数，收到 {measured_molar_mass}"
        )

    with localcontext() as ctx:
        # 内部精度远高于输出精度，避免摩尔比与反算百分比的中间舍入误差
        ctx.prec = 50

        total_pct = sum(percentages.values(), Decimal(0))
        if abs(total_pct - _HUNDRED) > tolerance:
            raise FormulaInferInputError(
                f"质量百分比之和为 {_round(total_pct, OUTPUT_PRECISION)}，"
                f"与 100 的偏差超过 tolerance={_round(tolerance, OUTPUT_PRECISION)}"
            )

        # ---- 质量百分比 -> 摩尔比 -> 归一 ----
        moles = {
            el: pct / ATOMIC_WEIGHTS[el] for el, pct in percentages.items()
        }
        n_min = min(moles.values())
        ratios = {el: n / n_min for el, n in moles.items()}

        # ---- 枚举 1～12 倍尺度的互质下标候选并择优 ----
        candidates = [
            c
            for scale in range(SCALE_MIN, SCALE_MAX + 1)
            if (c := _build_candidate(scale, ratios, percentages)) is not None
        ]
        best = _select_candidate(candidates, tolerance)

        subscripts: dict[str, int] = best["subscripts"]
        element_rows = [
            {
                "element": el,
                "subscript": subscripts[el],
                "given_percent": _round(percentages[el], OUTPUT_PRECISION),
                "theoretical_percent": _round(
                    best["theoretical"][el], OUTPUT_PRECISION
                ),
                "deviation": _round(best["deviations"][el], OUTPUT_PRECISION),
            }
            for el in percentages  # 按请求中的元素顺序回传
        ]

        result: dict[str, Any] = {
            "empirical_formula": best["hill_formula"],
            "subscripts": dict(subscripts),
            "formula_scale": best["scale"],
            "total_atoms": best["total_atoms"],
            "elements": element_rows,
            "max_abs_deviation": _round(
                best["max_abs_deviation"], OUTPUT_PRECISION
            ),
            "empirical_molar_mass": _round(
                best["molar_mass"], OUTPUT_PRECISION
            ),
            "measured_molar_mass": None,
            "molecular_multiplier": None,
            "molecular_formula": None,
            "molecular_molar_mass": None,
            "molar_mass_relative_deviation_percent": None,
        }

        # ---- 可选：实测摩尔质量 -> 分子式 ----
        if measured_molar_mass is not None:
            empirical_mass = best["molar_mass"]
            nearest = int(
                (measured_molar_mass / empirical_mass).quantize(
                    _ONE, rounding=ROUND_HALF_UP
                )
            )
            multiplier = max(1, nearest)  # 倍数必须为正整数
            molecular_mass = empirical_mass * multiplier
            relative_deviation = (
                abs(molecular_mass - measured_molar_mass)
                / measured_molar_mass
                * _HUNDRED
            )
            if relative_deviation > tolerance:
                raise NoMatchingFormulaError(
                    "实测摩尔质量与经验式最近正整数倍的相对偏差超过 "
                    "tolerance%，无法确定分子式",
                    reason="molar_mass_mismatch",
                    tolerance_percent=_round(tolerance, OUTPUT_PRECISION),
                    molecular_multiplier=multiplier,
                    theoretical_molar_mass=_round(
                        molecular_mass, OUTPUT_PRECISION
                    ),
                    relative_deviation_percent=_round(
                        relative_deviation, OUTPUT_PRECISION
                    ),
                )

            molecular_subscripts = {
                el: n * multiplier for el, n in subscripts.items()
            }
            result.update(
                {
                    "measured_molar_mass": _round(
                        measured_molar_mass, OUTPUT_PRECISION
                    ),
                    "molecular_multiplier": multiplier,
                    "molecular_formula": _format_formula(
                        molecular_subscripts, best["hill_order"]
                    ),
                    "molecular_molar_mass": _round(
                        molecular_mass, OUTPUT_PRECISION
                    ),
                    "molar_mass_relative_deviation_percent": _round(
                        relative_deviation, OUTPUT_PRECISION
                    ),
                }
            )

        return result
