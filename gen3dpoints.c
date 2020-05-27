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
//int plot()
{
  printf("x coordinates\n");
   int x[101];
   int counterx,countery,counterz; //counter variable
   FILE *fp;
   fp = fopen("coord.dat", "w+");
   if(fp==NULL)
   {
     printf("\nError:file cannot be opened");
     exit(1);
   }
   fprintf(fp, "This is testing for fprintf...\n");
   fprintf(fp,"x coordinates\n");
   //fputs("This is testing for fputs...\n", fp);
   //fclose(fp);
   for( counterx = 1; counterx <= 100; ++counterx )
   {
      x[counterx]= (rand()%100);
      printf("%5d", x[counterx]);
      fprintf(fp,"%5d",x[counterx]); //for the output in a file
      //for displaying value in new line
      if (counterx % 25 == 0)
      puts("");
      fputs("\n",fp);
   }
   
   fputs("\n",fp);
   //printf("%d\n",counterx);
   printf("%d\n",x[2]);

   printf("y coordinates\n");
   fprintf(fp,"y coordinates\n");
   int y[101];
   for (countery=1;countery<=100;++countery)
   {
     y[countery]=(rand()%100);
     printf("%5d",y[countery]);
     fprintf(fp,"%5d",y[countery]);
     if (countery%25==0)
     puts("");
     fputs("\n",fp);
   }
   fputs("\n",fp);
   printf("%d\n",y[2]);
   printf("8 days later, comment this later\n");
   printf("blah blah blah\n");

   printf("z coordinates\n");
   fprintf(fp,"z coordinates\n");
   int z[101];
   for (counterz=1;counterz<=100;++counterz)
   {
     z[counterz]=(rand()%100);
     printf("%5d",z[counterz]);
     fprintf(fp,"%5d",z[counterz]);
     if (counterz%25==0)
     puts("");
     fputs("\n",fp);
   }
   printf("%d\n",z[2]);
   fclose(fp);
   system("gnuplot -p file.gp");
   return 0;
   
}


