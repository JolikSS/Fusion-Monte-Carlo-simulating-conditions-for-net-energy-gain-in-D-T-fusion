#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// Та же функция, что на этапе 2. Позже вынесем общие функции в отдельный
double tunneling_probability(double E_keV, double EG_keV)
{
    return std::exp(-std::sqrt(EG_keV / E_keV));
}

int main()
{
    const double keV = 1.602176634e-16;    // 1 кэВ в джоулях
    const double m_D = 3.3435837724e-27;   // масса дейтрона, кг
    const double m_T = 5.0073567446e-27;   // масса тритона, кг
    const double kT_keV = 10.0;            // температура плазмы, кэВ
    const double EG_DT = 1182.0;           // энергия Гамова для D-T, кэВ
    const int N = 1'000'000;               // число смоделированных столкновений
    const unsigned SEED = 42;

    // Параметры гистограммы: энергии от 0 до 200 кэВ, корзины по 1 кэВ.
    // Столкновения с E > 200 кэВ отбрасываем: при kT = 10 кэВ их доля
    // порядка e^-20, на результат это не влияет.
    const double E_max = 200.0;
    const int n_bins = 200;
    const double bin_width = E_max / n_bins;

    // Приведённая масса: относительное движение двух частиц описывается
    // как движение одной частицы с массой m_r.
    const double m_r = m_D * m_T / (m_D + m_T);

    // Компоненты относительной скорости распределены нормально
    // с разбросом sqrt(kT / m_r), как на этапе 1, только масса другая.
    const double sigma = std::sqrt(kT_keV * keV / m_r);

    std::mt19937 rng(SEED);
    std::normal_distribution<double> normal(0.0, sigma);

    // collisions[i] - сколько столкновений попало в корзину i,
    // reactions[i]  - суммарный вес (вклад в число реакций) корзины i.
    // vector - массив, размер которого задаётся при создании.
    std::vector<double> collisions(n_bins, 0.0);
    std::vector<double> reactions(n_bins, 0.0);

    double sum_E = 0.0;   // для средней энергии всех столкновений
    double sum_w = 0.0;   // сумма всех весов (нормировка)
    double sum_wE = 0.0;  // для средней энергии, взвешенной по реакциям

    for (int i = 0; i < N; ++i) {
        double vx = normal(rng), vy = normal(rng), vz = normal(rng);
        double E = 0.5 * m_r * (vx * vx + vy * vy + vz * vz) / keV;

        // E = 0 теоретически возможно (дало бы деление на ноль в весе),
        // но вероятность этого ничтожна; пропускаем на всякий случай.
        if (E <= 0.0) continue;

        // Вклад столкновения в скорость реакции: sigma(E) * v ~ P(E) / sqrt(E).
        // Абсолютный множитель нам не нужен: сравниваем только форму распределений.
        double w = tunneling_probability(E, EG_DT) / std::sqrt(E);

        sum_E += E;
        sum_w += w;
        sum_wE += w * E;

        // Номер корзины: целая часть от E / ширина корзины.
        // static_cast<int> отбрасывает дробную часть (например, 12.7 -> 12).
        int bin = static_cast<int>(E / bin_width);
        if (bin < n_bins) {
            collisions[bin] += 1.0;
            reactions[bin] += w;
        }
    }

    std::ofstream out("gamow.csv");
    if (!out) {
        std::cerr << "не удалось открыть gamow.csv\n";
        return 1;
    }

    // Для каждой корзины записываем её центр и две доли:
    // какая часть всех столкновений и какая часть всех реакций приходится на неё.
    // Нормировка на N и sum_w делает обе колонки сопоставимыми (каждая в сумме ~1).
    out << "energy_keV,collisions_share,reactions_share\n";
    for (int b = 0; b < n_bins; ++b) {
        double E_center = (b + 0.5) * bin_width;
        out << E_center << ","
            << collisions[b] / N << ","
            << reactions[b] / sum_w << "\n";
    }

    // Теоретическое положение пика Гамова: E0 = (sqrt(EG) * kT / 2)^(2/3).
    // Это максимум произведения exp(-E/kT) * exp(-sqrt(EG/E)).
    double E0 = std::pow(std::sqrt(EG_DT) * kT_keV / 2.0, 2.0 / 3.0);

    std::cout << "Средняя энергия столкновения:           " << sum_E / N << " кэВ\n";
    std::cout << "Средняя энергия столкновения с реакцией: " << sum_wE / sum_w << " кэВ\n";
    std::cout << "Теоретический пик Гамова E0:            " << E0 << " кэВ\n";
    return 0;
}