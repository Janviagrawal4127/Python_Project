import numpy as np
import math

# Arithmetic Series
def arithmetic_series(first_term, common_diff, num_terms):
    """Calculate the sum of an arithmetic series."""
    series_sum = num_terms / 2 * (2 * first_term + (num_terms - 1) * common_diff)
    return series_sum

# Geometric Series
def geometric_series(first_term, common_ratio, num_terms):
    """Calculate the sum of a geometric series."""
    if common_ratio == 1:
        return first_term * num_terms
    series_sum = first_term * (1 - common_ratio**num_terms) / (1 - common_ratio)
    return series_sum

# Harmonic Series
def harmonic_series(num_terms):
    """Calculate the sum of a harmonic series."""
    series_sum = sum(1 / i for i in range(1, num_terms + 1))
    return series_sum

# Fourier Series Coefficients
def fourier_series_coeff(f, L, N):
    """Calculate the first N terms of the Fourier series coefficients for a function f over the interval [-L, L]."""
    dx = L * 1e-4
    x = np.arange(-L, L, dx)
    a0 = (1 / (2 * L)) * np.trapz([f(xi) for xi in x], x)
    an = [(1 / L) * np.trapz([f(xi) * np.cos(np.pi * n * xi / L) for xi in x], x) for n in range(1, N + 1)]
    bn = [(1 / L) * np.trapz([f(xi) * np.sin(np.pi * n * xi / L) for xi in x], x) for n in range(1, N + 1)]
    return a0, an, bn

# Taylor Series
def exp_taylor_series(x, terms=10):
    """Compute the Taylor series approximation of e^x using a specified number of terms."""
    result = sum(x**n / math.factorial(n) for n in range(terms))
    return result

# Binomial Series
def binomial_series(x, n):
    """Compute the binomial series expansion (1 + x)^n using the binomial theorem."""
    series_sum = 0
    for k in range(n + 1):
        binomial_coefficient = math.comb(n, k)  # Binomial coefficient: n choose k
        term = binomial_coefficient * (x ** k)
        series_sum += term
    return series_sum

# Alternating Series
def alternating_series_sum(start_term, num_terms):
    """Calculate the sum of an alternating series."""
    series_sum = 0
    sign = 1  # Start with positive sign
    for i in range(num_terms):
        term = sign * start_term  # Multiply the term with the current sign
        series_sum += term
        sign *= -1  # Change the sign for the next term
        start_term += 1  # Increment the starting term for the next iteration
    return series_sum

def main_menu():
    while True:
        print("\nMain Menu:")
        print("1. Arithmetic Series")
        print("2. Geometric Series")
        print("3. Harmonic Series")
        print("4. Fourier Series")
        print("5. Taylor Series")
        print("6. Binomial Series")
        print("7. Alternating Series")
        print("8. Exit")

        choice = input("Enter your choice: ")

        if choice == '1':
            first_term = float(input("Enter the first term: "))
            common_diff = float(input("Enter the common difference: "))
            num_terms = int(input("Enter the number of terms: "))
            print("Arithmetic Series Sum:", arithmetic_series(first_term, common_diff, num_terms))

        elif choice == '2':
            first_term = float(input("Enter the first term: "))
            common_ratio = float(input("Enter the common ratio: "))
            num_terms = int(input("Enter the number of terms: "))
            print("Geometric Series Sum:", geometric_series(first_term, common_ratio, num_terms))

        elif choice == '3':
            num_terms = int(input("Enter the number of terms: "))
            print("Harmonic Series Sum:", harmonic_series(num_terms))

        elif choice == '4':
            L = np.pi
            N = int(input("Enter the number of terms in the Fourier series: "))
            f = lambda x: np.abs(x)
            a0, an, bn = fourier_series_coeff(f, L, N)
            print("a0 =", a0)
            print("an =", an)
            print("bn =", bn)

        elif choice == '5':
            x_value = float(input("Enter the value of x to compute e^x: "))
            num_terms = int(input("Enter the number of terms in the Taylor series: "))
            approximation = exp_taylor_series(x_value, num_terms)
            exact_value = math.exp(x_value)
            print(f"Taylor Series approximation of e^{x_value} with {num_terms} terms:", approximation)
            print(f"Exact value of e^{x_value}:", exact_value)
            print(f"Error: {abs(exact_value - approximation)}")

        elif choice == '6':
            x_value = float(input("Enter the value of x in (1 + x)^n for binomial series: "))
            n_value = int(input("Enter the power n for the binomial series expansion: "))
            result = binomial_series(x_value, n_value)
            print(f"The binomial series sum for (1 + {x_value})^{n_value} is:", result)

        elif choice == '7':
            start_term = int(input("Enter the starting term of the alternating series: "))
            num_terms = int(input("Enter the number of terms in the alternating series: "))
            result = alternating_series_sum(start_term, num_terms)
            print("Alternating Series Sum:", result)

        elif choice == '8':
            print("Exiting the program... Thank you")
            break

        else:
            print("Invalid choice. Please enter a valid option.")

if __name__ == "__main__":
    main_menu()
