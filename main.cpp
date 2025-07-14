#include<iostream>
#include<stdlib.h>
#include<limits.h>

#include "EVRP.hpp"
#include "stats.hpp"
#include "genetic.hpp"

using namespace std;

/*initialiazes a run for your heuristic*/
void start_run(int r){
  srand(r); //random seed
  init_evals();
  init_current_best();
  cout << "Run: " << r << " with random seed " << r << endl;
}

/*gets an observation of the run for your heuristic*/
void end_run(int r){
  get_mean(r-1,get_current_best()); //from stats.h
  cout << "End of run " << r << " with best solution quality " << get_current_best() << " total evaluations: " << get_evals()  << endl;
  cout << " " << endl;
}

/****************************************************************/
/*                Main Function                                 */
/****************************************************************/
int main(int argc, char *argv[]) {

  if (argc < 2) {
    cout << "Usage: ./main <instance_file.evrp>" << endl;
    return 1;
  }

  int run;
  problem_instance = argv[1];       //pass the .evrp filename as an argument
  read_problem(problem_instance);   //Read EVRP from file from EVRP.h

  open_stats(); //open text files to store the best values from the 20 runs stats.h

  for(run = 1; run <= MAX_TRIALS; run++){
      start_run(run);

      // Substitui heurística por GA
      initialize_genetic(NUM_OF_CUSTOMERS);
      run_genetic();
      print_best_solution();
      end_run(run);  //store the best solution quality for each run
  }

  close_stats(); //close text files to calculate the mean result from the 20 runs stats.h

  //free memory
  free_stats();
  free_EVRP();

  return 0;
}
