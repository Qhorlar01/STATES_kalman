# STATES Project: Kalman Filter Implementation & Nonlinear Filters Comparison

This repository contains two projects focused on state estimation using Kalman filters and a comparison of nonlinear filtering techniques. The projects are part of the **Statistical Signal Processing and Estimation Theory (STATES)** course.

---

## Project 1: DC Motor State Estimation with Kalman Filter

### 📝 Description  
This project estimates the angular position $$(θ)$$ and angular velocity $$(Ω)$$ of a DC motor using **Stationary** and **Time-Varying Kalman Filters**. The system model includes a deterministic state-space representation and a stochastic linear model with sensor noise. Key tasks include:  
- Simulating input voltage (square wave).  
- Implementing Kalman filters for state estimation.  
- Tuning parameters (`q`, `r`, `P`) to analyze filter performance.  

### 🛠️ Implementation  
#### Key Functions:  
- `inputvoltage(D, A, $$Delta$$, $$Ts$$)`: Generates a zero-mean square wave input.  
- `simulate(u, G, T, Ts, L, $$x_1$$)`: Simulates the DC motor's deterministic model.  
- `stationary_kal(y, u, G, T, Ts, L, x1_θ, q)`: Implements the Stationary Kalman Filter.  
- `kal(y, u, G, T, Ts, L, x1_θ, p1_θ, q)`: Implements the Time-Varying Kalman Filter.  

#### Parameters:  
- Motor constants: `G = 50 rad·s⁻¹·V⁻¹`, `T = 0.02 s`  
- Sampling: `Ts = 0.001 s`, quantization levels `L`.  

### 📊 Results  
- **Effect of `q` (Process Noise Variance):**  
  - Small `q` (e.g., 0.01): Filter trusts the model, smooths sensor noise.  
  - Large `q` (e.g., 25): Filter relies on sensor data; velocity estimation becomes noisy.  
  - *Example Figure*:  
    ![Figure 6: Kalman Filter Performance for q=0.01 and q=25](figures/figure6.png)  

- **Model Robustness:**  
  - Filters tolerate parametric inaccuracies (e.g., `T_filter = 0.025 s` vs. `T_actual = 0.02 s`).  

---

## Project 2: Nonlinear Filters Comparison (EKF, UKF, Bootstrap Particle Filter)

### 📝 Description  
Compares the performance of **Extended Kalman Filter (EKF)**, **Unscented Kalman Filter (UKF)**, and **Bootstrap Particle Filter** for a nonlinear system:  
- **State Model**: `x[n+1] = 0.5x + 25x/(1+x²) + 8cos(1.2n) + v[n]`  
- **Observation Model**: `y[n] = x²/20 + w[n]`  
- Noise variances: `Q = 10`, `R = 1`.  

### 📊 Results  
- **Performance Ranking**:  
  1. **Bootstrap Particle Filter**: Best estimation accuracy.  
  2. **UKF**: Close performance with lower computational cost.  
  3. **EKF**: Diverges occasionally due to linearization errors.  
  - *Example Figure*:  
    ![Figure 16: State Estimation Comparison](figures/figure16.png)  

---

## 🚀 Getting Started  
### Requirements  
- MATLAB (tested on R2021a or later).  

### Repository Structure  
