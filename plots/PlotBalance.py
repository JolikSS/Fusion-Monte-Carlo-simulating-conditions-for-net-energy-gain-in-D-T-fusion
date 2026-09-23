import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/home/fuf/VisualCode/fusion_sim/src/balance.csv")

# Два графика рядом: слева баланс мощностей, справа критерий Лоусона.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# --- Баланс мощностей ---
ax1.plot(df["kT_keV"], df["P_fusion_MW"], "--", label="Вся мощность синтеза (17.6 МэВ)")
ax1.plot(df["kT_keV"], df["P_alpha_MW"], label="Нагрев альфа-частицами (3.5 МэВ)")
ax1.plot(df["kT_keV"], df["P_brems_MW"], label="Потери на излучение")

# Точка идеального зажигания: первая температура, где нагрев перекрывает излучение.
ignition = df[df["P_alpha_MW"] >= df["P_brems_MW"]]["kT_keV"].iloc[0]
ax1.axvline(ignition, color="gray", linestyle=":",
            label=f"Идеальное зажигание ≈ {ignition:.1f} кэВ")

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Температура kT, кэВ")
ax1.set_ylabel("Мощность, МВт/м³ (при n = 10²⁰ м⁻³)")
ax1.set_title("Нагрев против излучения")
ax1.legend(fontsize=8)
ax1.grid(True, which="both", alpha=0.4)

# --- Критерий Лоусона ---
# dropna() убирает температуры, где зажигание невозможно (там записано nan).
lawson = df.dropna(subset=["triple_product"])
ax2.plot(lawson["kT_keV"], lawson["triple_product"])

best = lawson.loc[lawson["triple_product"].idxmin()]
ax2.plot(best["kT_keV"], best["triple_product"], "o",
         label=f"Минимум: {best['kT_keV']:.1f} кэВ ≈ {best['kT_keV'] * 11.6:.0f} млн К")

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("Температура kT, кэВ")
ax2.set_ylabel("n·T·τ, кэВ·с/м³")
ax2.set_title("Критерий Лоусона: что нужно для зажигания")
ax2.legend()
ax2.grid(True, which="both", alpha=0.4)

plt.tight_layout()  # чтобы подписи двух графиков не налезали друг на друга
plt.savefig("balance.png", dpi=150)
plt.show()