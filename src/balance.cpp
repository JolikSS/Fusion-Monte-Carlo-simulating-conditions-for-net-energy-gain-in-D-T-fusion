#include <cmath>
#include <fstream>
#include <iostream>

// Скорость реакции <sigma v> для D-T по аппроксимации Bosch & Hale (1992),
// Nuclear Fusion 32, 611. Это подгонка экспериментальных данных, справедливая
// для 0.2-100 кэВ. Возвращает значение в м^3/с (сама формула даёт см^3/с).
double sigmav_DT(double T_keV)
{
    const double BG   = 34.3827;      // корень из энергии Гамова, кэВ^(1/2)
    const double mrc2 = 1124656.0;    // энергия покоя приведённой массы D-T, кэВ
    const double C1 = 1.17302e-9;
    const double C2 = 1.51361e-2;
    const double C3 = 7.51886e-2;
    const double C4 = 4.60643e-3;
    const double C5 = 1.35000e-2;
    const double C6 = -1.06750e-4;
    const double C7 = 1.36600e-5;

    double T = T_keV;
    // theta - "эффективная температура": поправка, которая учитывает
    // ядерный множитель (в том числе резонанс D-T), отсутствующий в нашей модели.
    double theta = T / (1.0 - T * (C2 + T * (C4 + T * C6))
                            / (1.0 + T * (C3 + T * (C5 + T * C7))));
    double xi = std::cbrt(BG * BG / (4.0 * theta));
    double sv_cm3 = C1 * theta * std::sqrt(xi / (mrc2 * T * T * T)) * std::exp(-3.0 * xi);
    return sv_cm3 * 1e-6;             // см^3 -> м^3
}

int main()
{
    const double keV_J   = 1.602176634e-16;   // 1 кэВ в джоулях
    const double E_alpha = 3.5e3 * keV_J;     // энергия альфа-частицы (3.5 МэВ) - греет плазму
    const double E_fus   = 17.6e3 * keV_J;    // полная энергия реакции (17.6 МэВ)
    const double C_B     = 5.35e-37;          // коэффициент тормозного излучения, Вт·м^3·кэВ^(-1/2)
    const double n       = 1e20;              // плотность ядер, м^-3 (типично для токамака)
    const double keV_to_MK = 11.6045;         // 1 кэВ = 11.6 миллиона кельвин

    std::ofstream out("balance.csv");
    if (!out) {
        std::cerr << "не удалось открыть balance.csv\n";
        return 1;
    }
    out << "kT_keV,P_fusion_MW,P_alpha_MW,P_brems_MW,triple_product\n";

    double T_ignition = -1.0;        // температура идеального зажигания (пока не найдена)
    double prev_ratio = 0.0, prev_T = 0.0;
    double best_T = 0.0, best_triple = 1e100;  // минимум тройного произведения

    for (double T = 1.0; T <= 100.0 + 1e-9; T += 0.05) {
        double sv = sigmav_DT(T);

        // Смесь 50/50: n_D = n_T = n/2, поэтому число реакций в секунду
        // в кубометре равно n_D * n_T * <sigma v> = n^2/4 * <sigma v>.
        double reactions = 0.25 * n * n * sv;

        double P_fusion = reactions * E_fus;              // вся выделяемая мощность, Вт/м^3
        double P_alpha  = reactions * E_alpha;            // мощность, остающаяся в плазме
        double P_brems  = C_B * n * n * std::sqrt(T);     // потери на излучение (n_e = n)

        // Идеальное зажигание: ищем, где P_alpha впервые превышает P_brems.
        // Точку пересечения уточняем линейной интерполяцией между шагами сетки.
        double ratio = P_alpha / P_brems;
        if (T_ignition < 0.0 && prev_T > 0.0 && prev_ratio < 1.0 && ratio >= 1.0) {
            T_ignition = prev_T + (1.0 - prev_ratio) * (T - prev_T) / (ratio - prev_ratio);
        }
        prev_ratio = ratio;
        prev_T = T;

        // Критерий Лоусона. Тепловая энергия плазмы в кубометре 3nkT
        // (ионы и электроны, по 3/2 kT на частицу). Она утекает за время tau_E,
        // и нагрев должен покрыть и излучение, и утечку:
        //     P_alpha - P_brems >= 3 n kT / tau_E.
        // Отсюда минимальное n*tau_E, а умножив на T, получаем тройное произведение.
        // Пока P_alpha < P_brems, зажигание невозможно ни при каком удержании: пишем nan.
        double triple = std::nan("");
        double net = 0.25 * sv * E_alpha - C_B * std::sqrt(T);   // (P_alpha - P_brems) / n^2
        if (net > 0.0) {
            double n_tau = 3.0 * T * keV_J / net;                  // с/м^3
            triple = n_tau * T;                                    // кэВ·с/м^3
            if (triple < best_triple) {
                best_triple = triple;
                best_T = T;
            }
        }

        // Мощности пишем в МВт/м^3 - так числа удобнее читать.
        out << T << "," << P_fusion / 1e6 << "," << P_alpha / 1e6 << ","
            << P_brems / 1e6 << "," << triple << "\n";
    }

    std::cout << "<sigma v> при 10 кэВ: " << sigmav_DT(10.0) << " м^3/с (справочно ~1.1e-22)\n\n";
    std::cout << "Идеальное зажигание: " << T_ignition << " кэВ = "
              << T_ignition * keV_to_MK << " млн К\n";
    std::cout << "Минимум тройного произведения: " << best_triple << " кэВ·с/м^3"
              << " при " << best_T << " кэВ = " << best_T * keV_to_MK << " млн К\n";
    return 0;
}