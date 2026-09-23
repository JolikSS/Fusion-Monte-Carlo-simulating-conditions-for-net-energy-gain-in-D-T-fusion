#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

double tunneling_probability(double E_keV, double EG_keV)
{
    return std::exp(-std::sqrt(EG_keV / E_keV));
}

// Скорость реакций <sigma v> методом Монте-Карло (с точностью до постоянного множителя).
// Генерируем N энергий столкновений при температуре kT и усредняем вес P(E)/sqrt(E).
// rng передаётся по ссылке (&), чтобы все вызовы продолжали одну и ту же
// последовательность случайных чисел, а не начинали её заново.
double rate_monte_carlo(double kT_keV, double EG_keV, int N, std::mt19937& rng)
{
    // Энергия частицы в распределении Максвелла имеет гамма-распределение
    // с параметром формы 3/2 и масштабом kT. Это то же распределение, что на этапе 1.
    std::gamma_distribution<double> maxwell_energy(1.5, kT_keV);

    double sum_w = 0.0;
    for (int i = 0; i < N; ++i) {
        double E = maxwell_energy(rng);
        if (E <= 0.0) continue;
        sum_w += tunneling_probability(E, EG_keV) / std::sqrt(E);
    }
    return sum_w / N;
}

// Аналитическое приближение для той же модели (метод перевала / гауссово
// приближение пика Гамова). Нужно для проверки моделирования.
double rate_theory(double kT_keV, double EG_keV)
{
    return std::pow(kT_keV, -2.0 / 3.0)
         * std::exp(-3.0 * std::cbrt(EG_keV / (4.0 * kT_keV)));
}

int main()
{
    const double EG_DT = 1182.0;
    const int N = 200'000;          // столкновений на каждую температуру
    const unsigned SEED = 42;
    const double kT_ref = 10.0;     // опорная температура: все значения делим на значение при 10 кэВ

    std::mt19937 rng(SEED);

    // Температуры от 1 до 100 кэВ, равномерно в логарифмическом масштабе:
    // kT = 10^(k/10), k = 0..20. Точка k = 10 даёт ровно 10 кэВ.
    // Логарифмическая сетка нужна потому, что на отрезке 1-10 кэВ происходит
    // всё самое интересное, а обычная сетка дала бы там слишком мало точек.
    std::vector<double> temps, mc, th;
    for (int k = 0; k <= 20; ++k) {
        double kT = std::pow(10.0, k / 10.0);
        temps.push_back(kT);
        mc.push_back(rate_monte_carlo(kT, EG_DT, N, rng));
        th.push_back(rate_theory(kT, EG_DT));
    }

    // Нормировочные значения при 10 кэВ (индекс 10 в массивах).
    double mc_ref = mc[10];
    double th_ref = th[10];

    std::ofstream out("rate.csv");
    if (!out) {
        std::cerr << "не удалось открыть rate.csv\n";
        return 1;
    }
    out << "kT_keV,rate_mc,rate_theory\n";

    std::cout << "  kT, кэВ    Монте-Карло      Теория\n";
    for (size_t i = 0; i < temps.size(); ++i) {
        double r_mc = mc[i] / mc_ref;
        double r_th = th[i] / th_ref;
        out << temps[i] << "," << r_mc << "," << r_th << "\n";
        std::cout << "  " << temps[i] << "    " << r_mc << "    " << r_th << "\n";
    }
    return 0;
}