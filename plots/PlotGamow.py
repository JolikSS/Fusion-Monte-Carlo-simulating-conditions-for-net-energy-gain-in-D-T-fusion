import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kT = 10.0
EG = 1182.0

df = pd.read_csv("/home/fuf/VisualCode/fusion_sim/src/gamow.csv")

# Две доли на одном графике: откуда берутся столкновения и откуда реакции.
# Обе колонки нормированы на 1, поэтому формы можно сравнивать напрямую.
plt.plot(df["energy_keV"], df["collisions_share"], label="Доля столкновений")
plt.plot(df["energy_keV"], df["reactions_share"], label="Доля реакций")

# Теоретическое положение пика Гамова для сверки с моделированием.
E0 = (np.sqrt(EG) * kT / 2) ** (2 / 3)
plt.axvline(E0, color="gray", linestyle="--", label=f"Пик Гамова, {E0:.1f} кэВ")

plt.xlim(0, 120)  # дальше обе кривые практически нулевые
plt.xlabel("Энергия столкновения, кэВ")
plt.ylabel("Доля на интервал 1 кэВ")
plt.title(f"Какие столкновения дают реакции D–T, kT = {kT} кэВ")
plt.legend()
plt.grid(True)
plt.savefig("gamow.png", dpi=150)
plt.show()