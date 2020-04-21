#include <random>

std::mt19937 eng; // object which produces random bits

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
eng.seed(seed); // seed the bit generator, replaces srand()

std::uniform_int_distribution<> dist(0, N-1); // encapsulates the correct method of turning
                                              // random bits into random numbers in the range
                                              // [0, N)
for(int i = 0; i < 6; i++) {
  int x = dist(eng); // the distribution internally used the engine like rand()
  int y = dist(eng);

  std::cout << "x:" << x << " y:" << y << "\n\n";
}