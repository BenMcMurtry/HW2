//  Created by Ben McMurtry on 15/10/2016.


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define RADIUS 0.310
#define DENSITY 1380
#define THICKNESS 0.00018
#define MIN 0
#define MAX 0.310
#define MAXEVAL 100
#define ABS_TOL 1e-15
#define REL_TOL 1e-15

typedef struct r{
    double root, uncertainty;
    int eval;
    struct r *next;
}Rootlist;

Rootlist *FirstOfList = NULL;
double (*InitialFunc)(double);


    
    
double rand_num (double min, double max) {
    double x, n, range;
    range = max = min;
    n = 2 * (double)(rand() / (double)RAND_MAX +1) - 1;
    x = min + (1 + n) * range / 2;
    return x;
}

void bisection(double(*func)(double), double *x1, double *x2, double *f1, double *f2) {
    double xNew, fNew;
    xNew = (*x1 + *x2) / 2;
    fNew = func(xNew);
    if (*f1 * fNew > 0) {
        *x1 = xNew;
        *f1 = fNew;
    } else {
        *x2 = xNew;
        *f2 = fNew;
    }
}

void secant(double(*func)(double), double *x1, double *x2, double *f1, double *f2) {
    double grad, intercept, xNew, fNew;
    grad = (*f2 - *f1) / (*x2 - *x1);
    intercept = (*f1 - grad * (*x1) + *f2 - grad * (*x2)) / 2;
    xNew = -intercept / grad;
    fNew = func(xNew);
    if (*f1 * fNew > 0) {
        *x1 = xNew;
        *f1 = fNew;
    } else {
        *x2 = xNew;
        *f2 = fNew;
    }
}

double defl_func(double x);

void sort_roots(void);

int solver(double (*func)(double), double guess1, double guess2, char method) {
    int found = 0, count = 0;
    double funcAns1 = func(guess1);
    double funcAns2 = func(guess2);
    //fn_def1 = fn(guess1) / deflating polynomial evaluated at guess1
    //fn_def2 = fn(guess2) / deflating polynomial evaluated at guess2
    if ( ! (funcAns1 * funcAns2 < 0.0) ) { /*if guess1 and guess2 have opposite sign*/
        switch (method) {
                case 'b' :
                    while ((fabs(guess1 - guess2) > (2 * ABS_TOL + REL_TOL * fabs(guess1 + guess2))) && (count < MAXEVAL)) {
                        bisection(func, &guess1, &guess2, &funcAns1, &funcAns2);
                        count++;
                }
                    FirstOfList = new_root(FirstOfList, guess1, guess2, count);
                    found++;
                    break;
                case 's' :
                    while ((fabs(guess1 - guess2) > (2 * ABS_TOL + REL_TOL * fabs(guess1 + guess2))) && (count < MAXEVAL)) {
                        secant(func, &guess1, &guess2, &funcAns1, &funcAns2);
                        count++;
                }
                    FirstOfList = new_root(FirstOfList, guess1, guess2, count);
                    found++;
                    break;
        }
    } else {
        double guess3, funcAns3;
        srand((unsigned int) time(NULL));
        do {
            guess3 = rand_num(guess2, guess1);
            funcAns3 = func(guess3);
            count++;
            if (count == MAXEVAL) {break;}
        } while(funcAns1 * funcAns3 >= 0);
        if (count < MAXEVAL) {
            found += solver(func, guess1, guess3, method);
            found += solver(func, guess3, guess2, method);
        }
    }
    return found;
}


int main() {
    double xMin, xMax;
    int found, foundcount;
    
    InitialFunc = &j0;
    xMin = 0;
    xMax = 100;
    
    found = solver(InitialFunc, xMin, xMax, 'b');
    foundcount = found;
    while(found != 0) {
        found = solver(&delf_func, xMin, xMax, 'b');
        foundcount = found;
    }
    
    if (foundcount ==0) {
        printf("\nNo roots were found between %g and %g", xMin, xMax);
    } else {
        sort_roots();
        printf("\n%d roots found\n", foundcount);
        for (Rootlist *p = FirstOfList; p != NULL; p = p->next) {
            printf("%g\t%g\t%d\n", p->root, p->uncertainty, p->eval);
        }
    }
    return 0;
}
