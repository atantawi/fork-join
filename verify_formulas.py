"""Verify the formulas in approximation-comparison.md"""
import sys
sys.path.insert(0, '/Users/tantawi/Projects/fork-join')

from forkjoin.analytical import (
    upper_bound_independent,
    lower_bound_bottleneck,
    nelson_tantawi,
    mean_response_time,
    mean_response_time_lh,
    _validate
)

def verify_homogeneous():
    """Verify homogeneous case formulas"""
    print("=" * 60)
    print("HOMOGENEOUS CASE VERIFICATION")
    print("=" * 60)

    mu = 1.0
    test_cases = [0.1, 0.3, 0.6, 0.9]

    for lam in test_cases:
        rho = lam / mu
        try:
            t_nelson = nelson_tantawi(lam, mu)
            t_ul = mean_response_time(lam, mu, mu)
            t_lh = mean_response_time_lh(lam, mu, mu)
            t_ub = upper_bound_independent(lam, mu, mu)
            t_bot = lower_bound_bottleneck(lam, mu, mu)

            print(f"\nlam={lam:.1f}, rho={rho:.2f}")
            print(f"  Nelson-Tantawi: {t_nelson:.6f}")
            print(f"  T_UL:            {t_ul:.6f}  (error: {(t_ul-t_nelson)/t_nelson*100:.4f}%)")
            print(f"  T_LH:            {t_lh:.6f}  (error: {(t_lh-t_nelson)/t_nelson*100:.4f}%)")
            print(f"  T_UB:            {t_ub:.6f}")
            print(f"  T_bot:           {t_bot:.6f}")

            # Verify T_UL is close to Nelson-Tantawi
            assert abs(t_ul - t_nelson) / t_nelson < 1e-10, "T_UL should exactly match Nelson-Tantawi"
            print("  ✓ T_UL = Nelson-Tantawi (exact)")

        except Exception as e:
            print(f"  Error: {e}")

def verify_heterogeneous():
    """Verify heterogeneous case formulas"""
    print("\n" + "=" * 60)
    print("HETEROGENEOUS CASE VERIFICATION")
    print("=" * 60)

    test_cases = [
        (1.0, 1.5, 0.3),
        (1.0, 1.5, 0.6),
        (1.0, 1.5, 0.9),
        (1.0, 2.0, 0.3),
        (1.0, 2.0, 0.6),
        (1.0, 2.0, 0.9),
        (1.0, 3.0, 0.6),
        (1.0, 5.0, 0.6),
    ]

    for mu1, mu2, lam in test_cases:
        try:
            rho1 = lam / mu1
            rho2 = lam / mu2
            t_ul = mean_response_time(lam, mu1, mu2)
            t_lh = mean_response_time_lh(lam, mu1, mu2)
            t_ub = upper_bound_independent(lam, mu1, mu2)
            t_bot = lower_bound_bottleneck(lam, mu1, mu2)

            print(f"\nμ1={mu1:.1f}, μ2={mu2:.1f}, λ={lam:.1f}, ρ1={rho1:.2f}, ρ2={rho2:.2f}")
            print(f"  T_UB:            {t_ub:.6f}")
            print(f"  T_bot:           {t_bot:.6f}")
            print(f"  T_UL:            {t_ul:.6f}")
            print(f"  T_LH:            {t_lh:.6f}")
            print(f"  T_UB >= T_UL >= T_bot: {t_ub >= t_ul >= t_bot}")

            # Verify bounds
            assert t_ub >= t_ul, "T_UB should be >= T_UL"
            assert t_ul >= t_bot, "T_UL should be >= T_bot"
            print("  ✓ Bounds satisfied")

        except Exception as e:
            print(f"  Error: {e}")

def verify_formulas():
    """Verify the simplified formulas match implementation"""
    print("\n" + "=" * 60)
    print("SIMPLIFIED FORMULA VERIFICATION")
    print("=" * 60)

    mu1 = 1.0
    mu2 = 2.0
    lam = 0.6

    # Implementation formula
    mu_min = min(mu1, mu2)
    mu_max = max(mu1, mu2)
    t0 = 1 / mu1 + 1 / mu2 - 1 / (mu1 + mu2)
    a0 = mu_min * t0
    alpha = mu_max / mu_min
    gamma = 1 + 0.375 * math.exp((1 - alpha) * 100)
    a1 = gamma - a0
    rho = lam / mu_min
    t_impl = (a0 + a1 * rho) / (mu_min - lam)

    # Simplified formula from document
    t_simplified = 1 / (mu_min - lam) + 1 / mu_max - 1 / (mu1 + mu2)

    print(f"\nμ1={mu1}, μ2={mu2}, λ={lam}")
    print(f"  Implementation:  {t_impl:.10f}")
    print(f"  Simplified:      {t_simplified:.10f}")
    print(f"  Difference:      {t_impl - t_simplified:.2e}")
    print(f"  Relative error:  {(t_impl - t_simplified)/t_simplified*100:.6f}%")

if __name__ == "__main__":
    import math
    verify_homogeneous()
    verify_heterogeneous()
    verify_formulas()