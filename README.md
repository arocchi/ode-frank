# ode-frank
ODEFrank is a small project that has the objective of bringing the ODE physics engine inside Simulink.

## Background
ODEFrank is a small piece of code I wrote during my coursework project at [Centro Piaggio](http://www.centropiaggio.unipi.it/) during my master thesis in 2012. The project consisted in implementing a bimanual peg in hole task with a modular humanoid robot (friendly called Frank) which was built using the VSA servo motors they developed. 

In the first phase of the project we needed a simple solution in order to implement  collision detection between the peg and the hole in three dimensions. This was achieved by using the [Open Dynamics Engine](http://www.ode.org/) physics engine. Since interfacing ODE within Simulink as a S-Function is a tedious practical matter, I decided to release it to the public. Even though now more modern method exist to simulate  robotics systems (gazebo, ROS interfaces for Matlab, yarp's Whole Body Interface (WBI)..), I decided to migrate the old google code project page to GitHub

The software has been tested under Windows and Linux, and is licensed under the Simplified BSD License.