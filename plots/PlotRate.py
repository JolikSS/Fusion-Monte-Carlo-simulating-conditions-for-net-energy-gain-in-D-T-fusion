import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/home/fuf/VisualCode/fusion_sim/src/rate.csv")

# Монте-Карло точками, теорию линией: так сразу видно, ложатся ли точки на кривую.
plt.plot(df["kT_keV"], df["rate_theory"], label="Аналитическое приближение")
plt.plot(df["kT_keV"], df["rate_mc"], "o", label="Монте-Карло")

# Обе оси логарифмические: температура меняется в 100 раз,
# а скорость реакций - на несколько порядков.
plt.xscale("log")
plt.yscale("log")

plt.xlabel("Температура kT, кэВ")
plt.ylabel("Скорость реакций относительно 10 кэВ")
plt.title("Зависимость скорости реакции D–T от температуры")
plt.legend()
plt.grid(True, which="both", alpha=0.4)  # which="both" - сетка и на промежуточных делениях
plt.savefig("../figures/rate.png", dpi=150)
plt.show()