"""经验式反推接口测试。

内置原子量（节选）：H 1.008、C 12.011、N 14.007、O 15.999、S 32.06、
Cl 35.45、K 39.0983、Fe 55.845
=> M(CH2O)=30.026, M(HO)=17.007, M(Fe2O3)=159.687
"""
from decimal import Decimal

import pytest
from fastapi.testclient import TestClient

from app.formula_infer import _build_candidate, _select_candidate
from app.main import app

client = TestClient(app)


def infer(elements, tolerance=None, molar_mass=None):
    payload = {"elements": elements}
    if tolerance is not None:
        payload["tolerance"] = tolerance
    if molar_mass is not None:
        payload["molar_mass"] = molar_mass
    return client.post("/api/v1/formula/infer", json=payload)


# ---------- 经验式反推 ----------

def test_fe2o3_requires_scale_2():
    # 69.94/30.06 的摩尔比约为 1:1.5，尺度 1（FeO）不达标，尺度 2 得 Fe2O3
    resp = infer({"Fe": "69.94", "O": "30.06"}, tolerance="0.1")
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["empirical_formula"] == "Fe2O3"
    assert data["subscripts"] == {"Fe": 2, "O": 3}
    assert data["formula_scale"] == 2
    assert data["total_atoms"] == 5
    assert Decimal(data["max_abs_deviation"]) <= Decimal("0.1")

    rows = {row["element"]: row for row in data["elements"]}
    assert rows["Fe"]["subscript"] == 2 and rows["O"]["subscript"] == 3
    assert Decimal(rows["Fe"]["theoretical_percent"]) == Decimal("69.943076")
    assert Decimal(rows["O"]["theoretical_percent"]) == Decimal("30.056924")
    # 偏差带符号：理论 - 实测
    assert Decimal(rows["Fe"]["deviation"]) == Decimal("0.003076")
    assert Decimal(rows["O"]["deviation"]) == Decimal("-0.003076")
    # 理论百分比合计 100
    total = sum(Decimal(row["theoretical_percent"]) for row in data["elements"])
    assert abs(total - 100) < Decimal("0.000002")

    # 未提供摩尔质量时分子式相关字段为 null
    assert data["molecular_formula"] is None
    assert data["molecular_multiplier"] is None
    assert data["molecular_molar_mass"] is None
    assert data["measured_molar_mass"] is None
    assert data["molar_mass_relative_deviation_percent"] is None
    assert data["empirical_molar_mass"] == "159.687000"


def test_glucose_empirical_and_molecular():
    resp = infer(
        {"H": "6.71", "C": "40.00", "O": "53.29"},
        tolerance="0.1",
        molar_mass="180.16",
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["empirical_formula"] == "CH2O"  # Hill：C、H 在前
    assert data["subscripts"] == {"C": 1, "H": 2, "O": 1}
    assert data["formula_scale"] == 1
    # 逐元素行保持请求顺序（H、C、O），与 Hill 化学式顺序解耦
    assert [row["element"] for row in data["elements"]] == ["H", "C", "O"]
    # M(CH2O)=30.026，最接近 180.16 的正整数倍为 6 倍 => C6H12O6
    assert data["molecular_multiplier"] == 6
    assert data["molecular_formula"] == "C6H12O6"
    assert data["molecular_molar_mass"] == "180.156000"
    assert Decimal(data["molar_mass_relative_deviation_percent"]) <= Decimal("0.1")
    assert data["measured_molar_mass"] == "180.160000"


def test_acetic_acid_molecular():
    resp = infer(
        {"C": "40.001", "H": "6.714", "O": "53.286"},
        molar_mass="60.05",
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["empirical_formula"] == "CH2O"
    assert data["molecular_multiplier"] == 2
    assert data["molecular_formula"] == "C2H4O2"
    assert Decimal(data["molar_mass_relative_deviation_percent"]) <= Decimal("0.1")


def test_h2o2_molecular():
    # 实测 34 g/mol：M(HO)=17.007，取 2 倍 => H2O2
    resp = infer({"H": "5.93", "O": "94.07"}, molar_mass="34")
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["empirical_formula"] == "HO"
    assert data["molecular_multiplier"] == 2
    assert data["molecular_formula"] == "H2O2"
    assert Decimal(data["molar_mass_relative_deviation_percent"]) <= Decimal("0.1")


def test_kcl_hill_ordering_without_carbon():
    # 不含碳时按元素符号字母序排列：Cl 在 K 前 => ClK
    resp = infer({"K": "52.45", "Cl": "47.55"})
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["empirical_formula"] == "ClK"
    assert data["subscripts"] == {"K": 1, "Cl": 1}
    assert [row["element"] for row in data["elements"]] == ["K", "Cl"]


def test_default_tolerance():
    resp = infer({"H": "6.71", "C": "40.00", "O": "53.29"})
    assert resp.status_code == 200, resp.text
    assert resp.json()["empirical_formula"] == "CH2O"


# ---------- 平局择优（直接覆盖纯函数选择规则） ----------

def _fake_candidate(formula, atoms, max_abs):
    return {
        "hill_formula": formula,
        "total_atoms": atoms,
        "max_abs_deviation": max_abs,
    }


def test_candidate_tie_break_fewer_atoms_then_hill_string():
    a = _fake_candidate("Fe2O3", 5, Decimal("0.05"))
    b = _fake_candidate("FeO", 2, Decimal("0.05"))
    assert _select_candidate([a, b], Decimal("0.1"))["hill_formula"] == "FeO"

    # 总原子数相同且最大偏差相同：取 Hill 字符串较小者（ClK < FeO）
    c = _fake_candidate("ClK", 2, Decimal("0.05"))
    assert _select_candidate([a, b, c], Decimal("0.1"))["hill_formula"] == "ClK"

    # 最大偏差不同：偏差最小者胜，即使原子数更多
    d = _fake_candidate("Fe2O3", 5, Decimal("0.01"))
    assert _select_candidate([b, c, d], Decimal("0.1"))["hill_formula"] == "Fe2O3"


def test_non_coprime_scales_are_dropped():
    # 摩尔比恰好 2:1：仅尺度 1 互质；尺度 2（4:2，gcd=2）必须被丢弃
    c1 = _build_candidate(1, {"H": Decimal(2), "O": Decimal(1)}, {"H": Decimal(1), "O": Decimal(99)})
    c2 = _build_candidate(2, {"H": Decimal(2), "O": Decimal(1)}, {"H": Decimal(1), "O": Decimal(99)})
    assert c1 is not None and c1["subscripts"] == {"H": 2, "O": 1}
    assert c2 is None


# ---------- 422 ----------

def test_no_empirical_candidate_returns_422():
    # 摩尔比接近黄金分割：1～12 尺度内最佳近似 C13N8 的最大偏差仍 > 0.1
    resp = infer({"C": "58.1145", "N": "41.8855"}, tolerance="0.1")
    assert resp.status_code == 422, resp.text
    detail = resp.json()["detail"]
    assert detail["reason"] == "no_empirical_candidate"
    assert detail["best_formula"] == "C13N8"
    assert Decimal(detail["best_max_abs_deviation"]) > Decimal("0.1")


def test_molar_mass_relative_deviation_too_large_returns_422():
    # CH2O 经验式 30.026，实测 145：最近倍数 5 => 150.13，相对偏差 3.54% > 0.1%
    resp = infer(
        {"H": "6.71", "C": "40.00", "O": "53.29"},
        tolerance="0.1",
        molar_mass="145",
    )
    assert resp.status_code == 422, resp.text
    detail = resp.json()["detail"]
    assert detail["reason"] == "molar_mass_mismatch"
    assert detail["molecular_multiplier"] == 5
    assert Decimal(detail["relative_deviation_percent"]) > Decimal("0.1")


# ---------- 400 / 模型校验 422 ----------

def test_non_positive_percent_returns_400():
    resp = infer({"H": "0", "O": "100"})
    assert resp.status_code == 400
    assert "正数" in resp.json()["detail"]


def test_negative_percent_returns_400():
    resp = infer({"H": "-1", "O": "101"})
    assert resp.status_code == 400


def test_sum_deviation_beyond_tolerance_returns_400():
    resp = infer({"H": "40", "O": "50"}, tolerance="0.1")
    assert resp.status_code == 400
    assert "100" in resp.json()["detail"]


def test_sum_deviation_within_tolerance_passes():
    # 合计 100.05，容差 0.1，边界内允许
    resp = infer({"H": "11.24", "O": "88.81"}, tolerance="0.1")
    assert resp.status_code == 200, resp.text
    assert resp.json()["empirical_formula"] == "H2O"


def test_unknown_element_returns_400():
    resp = infer({"Zn": "50", "O": "50"})
    assert resp.status_code == 400
    assert "Zn" in resp.json()["detail"]
    assert "原子量" in resp.json()["detail"]


def test_empty_elements_returns_400():
    resp = client.post("/api/v1/formula/infer", json={"elements": {}})
    assert resp.status_code == 400
    assert "至少" in resp.json()["detail"]


def test_invalid_tolerance_returns_422():
    resp = infer({"H": "11.2", "O": "88.8"}, tolerance="0")
    assert resp.status_code == 422


def test_non_positive_molar_mass_returns_400():
    resp = infer({"H": "11.2", "O": "88.8"}, molar_mass="0")
    assert resp.status_code == 400
    assert "molar_mass" in resp.json()["detail"]


@pytest.mark.parametrize(
    "payload",
    [
        {"tolerance": "0.1"},  # 缺 elements
        {"elements": {"H": "11.2", "O": "88.8"}, "molar_mass": "abc"},
    ],
)
def test_malformed_payload_returns_422(payload):
    resp = client.post("/api/v1/formula/infer", json=payload)
    assert resp.status_code == 422
