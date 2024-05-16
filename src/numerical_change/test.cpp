#include <iostream>
#include <cmath>

enum Sex {
    Male = 0,
    Female = 1
};

int main()
{
    for (int dim_i = 0; dim_i < 4; ++dim_i)
    {
        bool envt = static_cast<int>(std::floor(dim_i/2.0));
        Sex s1 = static_cast<Sex>(dim_i % 2);

        std::cout << envt << " " << s1 << std::endl;
    }

}
