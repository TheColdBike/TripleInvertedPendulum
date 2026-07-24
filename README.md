# TripleInvertedPendulum
A simulated Triple inverted pendulum controls project with noise and Kalman Filtering. 

# Controls Methods Used
An LQR is used on a linearized model from the upright position of the triple pendulum. This linearized model was found from the full dynamics of the triple pendulum on a cart, which were derived using Lagrangian Mechanics. 

# State Estimation Method Used
An Unscented Kalman Filter is used to estimate the state of the triple inverted pendulum. It is assumed that sensors for the angle of each pendulum as well as the position of the cart are present. Noise is injected into the sensors to ensure the true state is not known to the estimator. 
