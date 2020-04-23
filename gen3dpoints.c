/*#include <random>

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
*/

/*#include <random>
#include <iostream>
#include <vector>
#include <string>

int main() {
    const int N = 10;
    std::vector<std::string> square(N, std::string(N, '*'));
    
    for (auto &&row : square) {
        std::cout << row << '\n';
    }
    
    std::mt19937 eng;
    eng.seed(std::random_device{}());

    std::uniform_int_distribution<> dist(0, N-1);
    
    for(int i = 0; i < 6; i++) {
      int x = dist(eng);
      int y = dist(eng);
      
      std::cout << "x:" << x << " y:" << y << "\n\n";
      
      square[y][x] = ' ';
    }
    
    for (auto &&row : square) {
        std::cout << row << '\n';
    }


}
*/
//test
#include<stdio.h>
#include<stdlib.h>

int main()
{
    printf("x coordinates\n");
   int counterx,countery,counterz; //counter variable
   
   //displaying random numbers between 1 and 6
   for( counterx = 1,countery=1,counterz=1; counterx <= 100,countery<=100,counterz<=100; ++counterx,++countery,++counterz )
   {
      printf("x");
      printf("%5d", ( rand() % 100 ));
      printf("y");
      printf("%5d",(rand()%100));
      printf("z");
      printf("%5d",(rand()%100));

      //for displaying value in new line
      if (counterx,countery,counterz % 25 == 0)
      puts("");
   }
   printf("y coordinates\n");
   
   return 0;
   
}


