#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct EQUATION
{
    int y_coefficient;
    int y_prime_coefficient;
    int gx_term_count;
    float *gx_coefficients;
    int *gx_powers;
} EQUATION;

void printAbsoluteError(float trueValue, float approximateValue);

void printEquation(EQUATION *equation);

float calculateDerivative(float x, float y, EQUATION *equation);

float rungeKutta(float x0, float y0, float x_target, float h, EQUATION *equation);

float calculate_gx(EQUATION *equation, float t);

int main()
{
    int i;
    float x0, y0, h, x_target;
    float trueValue;

    EQUATION *equation = (EQUATION *)malloc(sizeof(EQUATION));
    if (equation == NULL) {
        perror("Memory allocation failed for equation");
        return 1;
    }

    printf("Runge-Kutta method for first-order differential equation\n");
    printf("a * y' + b * y = g(x)\n\n");

    printf("Enter the coefficient of y': ");
    scanf("%d", &equation->y_prime_coefficient);
    printf("Enter the coefficient of y: ");
    scanf("%d", &equation->y_coefficient);
    printf("Enter the number of terms in g(x): ");
    scanf("%d", &equation->gx_term_count);

    equation->gx_coefficients = (float *)malloc(equation->gx_term_count * sizeof(float));
    equation->gx_powers = (int *)malloc(equation->gx_term_count * sizeof(int));

    if (equation->gx_coefficients == NULL || equation->gx_powers == NULL) {
         perror("Memory allocation failed for g(x) terms");
         free(equation);
         return 1;
    }


    for (i = 0; i < equation->gx_term_count; i++)
    {
        printf("Enter the coefficient for the %d. term of g(x): ", i + 1);
        scanf("%f", &equation->gx_coefficients[i]);
        printf("Enter the power for the %d. term of g(x): ", i + 1);
        scanf("%d", &equation->gx_powers[i]);
    }

    printf("Entered equation: ");
    printEquation(equation);

    printf("Enter initial values for the solution:\n");
    printf("x0: ");
    scanf("%f", &x0);
    printf("y0: ");
    scanf("%f", &y0);

    printf("Enter the step size (h): ");
    scanf("%f", &h);

    printf("Enter the x value to find the solution at: ");
    scanf("%f", &x_target);

    printf("Enter the true value of the function at the target point: ");
    scanf("%f", &trueValue);

    printf("\n\n");
    float result = rungeKutta(x0, y0, x_target, h, equation);

    printf("y(%f) = %f\n", x_target, result);
    printf("Approximate value of the equation at point %f: %f\n", x_target, result);

    printAbsoluteError(trueValue, result);

    free(equation->gx_coefficients);
    free(equation->gx_powers);
    free(equation);

    return 0;
}

void printEquation(EQUATION *equation)
{
    int i;
    printf("%d*y' + %d*y = ", equation->y_prime_coefficient, equation->y_coefficient);
    for (i = 0; i < equation->gx_term_count; i++)
    {
        printf("( %lf * x^%d )%s", equation->gx_coefficients[i], equation->gx_powers[i], (i < equation->gx_term_count - 1) ? " + " : "");
    }
    printf("\n");
}

float calculateDerivative(float x, float y, EQUATION *equation)
{
    if (equation->y_prime_coefficient == 0) {
        fprintf(stderr, "Error: Coefficient of y' is zero, cannot calculate derivative.\n");
        exit(EXIT_FAILURE);
    }
    return (calculate_gx(equation, x) - equation->y_coefficient * y) / equation->y_prime_coefficient;
}

float rungeKutta(float x0, float y0, float x_target, float h, EQUATION *equation)
{
    int i;
    float k1, k2, k3, k4;
    float y = y0;
    float current_x = x0;

    if (x_target < x0)
    {
        h = -h;
    }

    int n = (int)((x_target - x0) / h);

    printf("ITERATIONS\n");
    for (i = 0; i < n; i++)
    {
        printf("f(%lf)=%lf\n", current_x, y);
        k1 = h * calculateDerivative(current_x, y, equation);
        k2 = h * calculateDerivative(current_x + 0.5 * h, y + 0.5 * k1, equation);
        k3 = h * calculateDerivative(current_x + 0.5 * h, y + 0.5 * k2, equation);
        k4 = h * calculateDerivative(current_x + h, y + k3, equation);

        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        current_x = current_x + h;
    }

    return y;
}

float calculate_gx(EQUATION *equation, float t)
{
    int i;
    float result = 0;
    for (i = 0; i < equation->gx_term_count; i++)
    {
        result += equation->gx_coefficients[i] * pow(t, equation->gx_powers[i]);
    }
    return result;
}

void printAbsoluteError(float trueValue, float approximateValue)
{
    printf("Absolute error: %f\n", fabs(trueValue - approximateValue));
}
