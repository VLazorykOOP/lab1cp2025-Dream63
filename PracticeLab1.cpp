#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

struct Entry {
    double x;
    double t;
    double u;
};

void printTable(const std::vector<Entry>& data) {
    std::cout << std::setw(10) << "X"
        << std::setw(10) << "T"
        << std::setw(10) << "U" << '\n';

    std::cout << std::string(30, '-') << '\n';

    for (const auto& entry : data) {
        std::cout << std::setw(10) << std::fixed << std::setprecision(3) << entry.x
            << std::setw(10) << entry.t
            << std::setw(10) << entry.u << '\n';
    }
}

std::string errorMessage1 = "Unable to open file: ",
errorMessage2 = "z * z + x * y < 0. Recalc Rrz with algorithm2.",
errorMessage3 = "x * x + z * y < 0. Recalc Rrz with algorithm3.",
errorMessage4 = "x * x + z * y < 0. Recalc Krn with algorithm4.";


std::vector<Entry> loadData(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::string msg = errorMessage1 + filename;
        throw std::runtime_error(msg);
    }

    std::vector<Entry> data;
    std::string line;

    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string xs, ts, us;
        std::getline(ss, xs, ',');
        std::getline(ss, ts, ',');
        std::getline(ss, us, ',');

        Entry entry;
        entry.x = std::stod(xs);
        entry.t = std::stod(ts);
        entry.u = std::stod(us);
        data.push_back(entry);
    }

    return data;
}
std::vector<Entry> data;

bool isBetween(double a, double b, double x) {
    return (a <= x && x <= b) || (b <= x && x <= a);
}
double getU(double x) {
    if (abs(x) <= 1) data = loadData("dat_X_1_1.dat");
    else if (x < -1) {
        x = -1 / x;
        data = loadData("dat_X_1_00.dat");
    }
    else if (x > 1) {
        x = -1 / x;
        data = loadData("dat_X_00_1.dat");
    }

    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i].x == x) {
            return data[i].u;
        }
        if (isBetween(1, data.size() - 2, i))
            if (isBetween(data[i - 1].x, data[i].x, x))
                return data[i - 1].u + (data[i].u - data[i].u) * (x - data[i - 1].x) / (data[i].x - x);
    }
    x = 0;
    throw std::out_of_range("x is out of data range");
}
double getT(double x) {
    if (abs(x) <= 1) data = loadData("dat_X_1_1.dat");
    else if (x < -1) {
        x = -1 / x;
        data = loadData("dat_X_1_00.dat");
    }
    else if (x > 1) {
        x = -1 / x;
        data = loadData("dat_X_00_1.dat");
    }

    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i].x == x) {
            return data[i].t;
        }
        if (isBetween(1, data.size() - 2, i))
            if (isBetween(data[i - 1].x, data[i].x, x))
                return data[i-1].t + (data[i].t - data[i].t) * (x - data[i-1].x) / (data[i].x - x);
    }

    throw std::out_of_range("x is out of data range");
}

double a1_Srz(double x, double y, double z) {
    if (x > y) return getT(x) + getU(z) - getT(y);
    return getU(y) + getU(y) - getU(z);
}
double a1_Srs1(double x, double y, double z) {
    if (z > y) {
        if (z * z + x * y < 0) throw std::runtime_error(errorMessage2);
        return a1_Srz(x, y, z) + y * log(z * z + x * y);
    }
    if (x * x + z * y < 0) throw std::runtime_error(errorMessage3);
    return y + a1_Srz(x, y, z) * sqrt(x * x + z * y);
}
double a1_Srs(double x, double y, double z) {
    if (z > y) {
        if (z * z + x * y <= 1) throw std::runtime_error(errorMessage2);
        return a1_Srz(x, y, z) + y * sqrt(z * z + x * y);
    }
    if (x * x + z * y <= 1) throw std::runtime_error(errorMessage4);
    return y + a1_Srz(z, x, y) * sqrt(x * x + z * y);
}
double a1_Qrz(double x, double y) {
    if (abs(x) < 1) return x * a1_Srs(x, y, x);
    return y * a1_Srs1(y, x, y);
}
double a1_Rrz(double x, double y, double z) {
    if (x > y) return x * z * a1_Qrz(y, z) - x;
    return y * x * a1_Qrz(x, y) + y;
}
double a1_Krn(double x, double y, double z) {
    return 73.1389 * a1_Rrz(x, y, y) + 14.838 * a1_Rrz(x - y, z, y);
}
double a1_fun(double x, double y, double z) {
    return x * a1_Krn(x, y, z) + y * a1_Krn(x, z, y) - z * a1_Krn(x, z, y);
}

double a2_Srs(double x, double y, double z) {
    if (z > y) return a1_Srz(x, y, z) + 1.44 * y * z;
    return y + 1.44 * a1_Srz(z, x, y);
}
double a2_Qrz(double x, double y) {

    if (abs(x) < 1) return x * a2_Srs(x, y, x);
    return y * a2_Srs(y, x, y);
}
double a2_Rrz(double x, double y, double z) {

    if (x > y) return x * y * a2_Qrz(y, z);
    return x * z * a2_Qrz(x, y);
}
double a2_Krn(double x, double y, double z) {

    return 73.1389 * a2_Rrz(x, y, y) + 14.838 * a2_Rrz(x - y, z, y);
}
double a2_fun(double x, double y, double z) {

    return x * a2_Krn(x, y, z) + y * a2_Krn(x, z, y) - z * a2_Krn(x, z, y);
}

double a3_Srs(double x, double y, double z) {
    if (z > y) return a1_Srz(x, y, z) + y * x;
    return y * z + a1_Srz(z, x, y);
}
double a3_Qrz(double x, double y) {
    if (abs(x) < 1) return x * a3_Srs(x, y, x);
    return y * a3_Srs(y, x, y);
}
double a3_Rrz(double x, double y, double z) {
    if (x > y) return x * y * a3_Qrz(y, z);
    return y * z * a2_Qrz(x, y);
}
double a3_Krn(double x, double y, double z) {
    return 73.1389 * a3_Rrz(x, y, y) + 14.838 * a3_Rrz(x - y, z, y);
}
double a3_fun(double x, double y, double z) {
    return x * a3_Krn(x, y, z) + y * a3_Krn(x, z, y) - z * a3_Krn(x, z, y);
}

double a4_Srs(double x, double y, double z) {
    if (z > y) return a1_Srz(x, y, z) + y * x;
    return y * z + a1_Srz(z, x, y);
}
double a4_Qrz(double x, double y) {
    if (abs(x) < 1) return x * a4_Srs(x, y, x);
    return x * a4_Srs(y, x, y);
}
double a4_Rrz(double x, double y, double z) {
    if (x > y) return y * a4_Qrz(y, z);
    return z * a4_Qrz(x, y);
}
double a4_Krn(double x, double y, double z) {
    return 83.1389 * a4_Rrz(x, y, z) + 4.838 * a4_Rrz(x, z, y);
}
double a4_fun(double x, double y, double z) {
    return x * a4_Krn(x, y, z) + y * a4_Krn(x, z, y) - z * a4_Krn(x, z, y);
}

double a5_fun(double x, double y, double z) {
    return 4.349 * x * z + 23.23 * y - 2.348 * x * y * z;
}

int main()
{
    std::string filename;
    double x, y, z;
    int alg = 1;
    std::cout << "Enter x: ";
    std::cin >> x;
    std::cout << "Enter y: ";
    std::cin >> y;
    std::cout << "Enter z: ";
    std::cin >> z;

    double value = 0;


    // use algorithm depending on some values
    algs:try {
        switch (alg)
        {
        case 1:
            value = a1_fun(x, y, z);
            break;
        case 2:
            value = a2_fun(x, y, z);
            break;
        case 3:
            value = a3_fun(x, y, z);
            break;
        case 4:
            value = a4_fun(x, y, z);
            break;
        case 5:
            value = a5_fun(x, y, z);
            break;
        default:
            break;
        }
    }
    catch (const std::runtime_error& ex)
    {
        std::cout << ex.what();
        
        if (((std::string)ex.what()).find(errorMessage1)!= std::string::npos)
            alg = 5;

        if (ex.what() == errorMessage2) 
            alg = 2;
        
        if (ex.what() == errorMessage3)
            alg = 3;

        if (ex.what() == errorMessage4)
            alg = 4;

        goto algs;
    }

    std::cout << std::endl << "Fun value: " << value;
}