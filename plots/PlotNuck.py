import numpy as np # для вычисления теоретической кривой
import pandas as pd # для чтения CSV-файла
import matplotlib.pyplot as plt # для построения графиков

kT = 10.0 # температура в кэВ, та же, что в программе на C++

E = pd.read_csv("/home/fuf/VisualCode/fusion_sim/src/maxwell.csv")["energy_keV"] # читаем энергии частиц

# гистограмма энергий из моделирования
plt.hist(E, bins=200, density=True, alpha=0.6, label="Моделирование")

# теоретическое распределение Максвелла по энергиям
x = np.linspace(0, E.max(), 500)
f = 2 * np.sqrt(x / np.pi) * kT ** -1.5 * np.exp(-x / kT)
plt.plot(x, f, "r", label="Теория")

plt.xlabel("Энергия частицы, кэВ") # подпись оси x
plt.ylabel("Плотность вероятности") # подпись оси y
plt.title(f"Распределение Максвелла, kT = {kT} кэВ")
plt.legend()
plt.savefig("maxwell.png", dpi=150) # сохраняем картинку
plt.show() # показываем график