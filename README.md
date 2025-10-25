# Lotka-Volterra-Competition-MATLAB
# 🦋 Population Dynamics of Two Competing Species

This MATLAB project models the population dynamics of two interacting species using the classic **Lotka-Volterra competition model**. The system of coupled ordinary differential equations (ODEs) is solved numerically using a custom-built **4th Order Runge-Kutta (RK4) method**.

The script generates phase portraits and time-evolution plots for four distinct ecological scenarios, providing a clear visual understanding of how species interactions and initial conditions determine their long-term survival.



---

## 🔬 The Mathematical Model

The simulation is based on the following system of non-linear ODEs that describe two competing species:

$$
\frac{dN_1}{dt} = r_1 N_1 \left(1 - \frac{N_1}{K_1}\right) - b_1 N_1 N_2
$$
$$
\frac{dN_2}{dt} = r_2 N_2 \left(1 - \frac{N_2}{K_2}\right) - b_2 N_1 N_2
$$

Where:
* $N_1, N_2$ are the populations of species 1 and 2.
* $r_1, r_2$ are the intrinsic growth rates.
* $K_1, K_2$ are the environmental carrying capacities.
* $b_1, b_2$ are the competition coefficients.

For a more general analysis, the code non-dimensionalizes these equations into the following form:

$$
\frac{dx}{d\tau} = x(1 - x) - \beta_1 xy
$$
$$
\frac{dy}{d\tau} = \alpha y \left(1 - \frac{y}{\kappa}\right) - \beta_2 xy
$$

This allows us to study the system's behavior based on key parameter ratios ($\alpha, \beta_1, \beta_2, \kappa$) rather than specific dimensional values.

---

## ✨ Features & Scenarios Simulated

The script automatically runs and visualizes four fundamental outcomes of species competition:

1.  **Stable Coexistence:** Both species survive and maintain stable populations over time.
2.  **Competitive Exclusion (Species 1 Wins):** Species 1 outcompetes and drives Species 2 to extinction.
3.  **Competitive Exclusion (Species 2 Wins):** Species 2 outcompetes and drives Species 1 to extinction.
4.  **Bistability:** The outcome (which species wins) depends entirely on the initial population sizes.

For each case, the script generates:
* A **Phase Portrait** showing multiple population trajectories and the system's nullclines.
* A **Time Series Plot** showing how the populations evolve over time from a specific starting point.

---

## 🚀 How to Use

1.  **Prerequisites:** Make sure you have **MATLAB** installed.
2.  **Download:** Download the `population_dynamics.m` file from this repository.
3.  **Run:** Open the file in MATLAB and click the **Run** button in the editor tab.

The script will execute and automatically generate four separate figure windows, one for each of the scenarios described above.

---

## ⚙️ Code Structure

The code is organized into four main functions:

* **`population_dynamics`**: The main driver function. It defines the model parameters and calls the plotting function for each of the four cases.
* **`plotPhasePortrait`**: Takes a set of non-dimensional parameters and generates the phase portrait and time series plots for that case.
* **`rk4_solver`**: The core numerical engine. It solves the system of ODEs using the 4th Order Runge-Kutta algorithm.
* **`f`**: Defines the non-dimensional system of ODEs to be solved.

Feel free to modify the parameters in the `population_dynamics` function to explore different system behaviors!
