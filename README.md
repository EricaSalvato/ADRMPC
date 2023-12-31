# ADRMPC

The folder contains the code used to obtain the results presented  in the paper:

Salvato E. , Fenu G. , Pellegrino F. A., Parisini T. "*An Active Disturbance Rejection Model Predictive Controller for Constrained Over-actuated Systems,* " submitted to the European Control Conference 2024, ECC24




The MATLAB code allows to apply ADRMPC to a simplified modelization of synchrotron, i.e. 
assuming that the dynamic model of the facility is defined as follows:

$$
\left\{ \begin{array}{l}
\dot{x}(t) = A x(t) -Au(t) \\
y(t) =  R x(t) + d(t) 
\end{array} \right.
$$


where $x \in \mathbb{R}^n$ is the state corresponding to the correction channels ($n=m$), $d$ represents both measurement noise and the disturbances, $A=\text{diag}(-\lambda_1,-\lambda_2,\ldots,-\lambda_n)$ describes the low-pass dynamics of the correction channels, where each $\lambda_i$, $i=1,2,\ldots ,n$, depends on the cut-off frequency value of the corresponding corrector magnet and $R$ assumes the following value: 

 $$ R=\begin{bmatrix}
     -0.307476818561554	&-0.632372349500656	&-0.454434603452683	&0.476383119821548	&0.406481325626373	&0.588919669389725	&0.463898330926895	&-0.689748972654343	&-0.508117973804474	&-0.240872204303741	&-0.188872218132019	&-0.653477013111115	&0.177061855792999	&-0.691902041435242	&0.645536184310913\\
 0.514381397515535	&0.385068040341139	&0.210438389331102	&-0.118657415732741	&-0.166663266718388	&-0.440627802163363	&-0.395499672740698	&0.443446002900600	&0.514566581696272	&0.324061121791601	&-0.226019881665707	&0.289400257170200	&0.180434379726648	&0.477587264031172	&-0.274364631623030\\
 -0.640640854835510	&-0.291974544525147	&-0.0731462240219116	&0.380523800849915	&0.359310507774353	&0.610096454620361	&0.504031479358673	&-0.694168806076050	&-0.598398447036743	&-0.324738323688507	&-0.0396096706390381	&-0.610961318016052	&0.0296482443809509	&-0.718413591384888	&0.570723414421082\\
 -0.417208224534988	&-0.353795886039734	&-0.191925317049027	&0.322388410568237	&0.212885886430740	&0.308444201946259	&0.288768857717514	&-0.326716154813767	&-0.427525341510773	&-0.285479873418808	&0.221302956342697	&-0.214967131614685	&-0.239037722349167	&-0.390224605798721	&0.129879415035248\\
 0.323294624686241	&0.723417401313782	&0.537412315607071	&-0.616078749299049	&-0.539144873619080	&-0.764655470848084	&-0.626091882586479	&0.867292284965515	&0.647285729646683	&0.422865971922874	&-0.879719406366348	&-0.355139225721359	&0.862076133489609	&0.0762896239757538	&0.521742850542069\\
 0.173676796257496	&0.438016541302204	&0.333082601428032	&-0.367206595838070	&-0.328615680336952	&-0.488302260637283	&-0.404852479696274	&0.546820685267448	&0.430314615368843	&0.251107625663281	&-0.547757726162672	&-0.244370102882385	&0.542351678013802	&0.0170819461345673	&0.352637432515621\\
 -0.490697994828224	&-0.682480707764626	&-0.458479151129723	&0.575155615806580	&0.454199537634850	&0.485072657465935	&0.367260202765465	&-0.635165944695473	&-0.320197716355324	&-0.114135965704918	&0.663314983248711	&0.0361488014459610	&-0.707611069083214	&-0.345796719193459	&-0.240701138973236
 \end{bmatrix}$$ 

* * *
## Results



![ADRMPC_outputsVSrefs_inputs.svg](_resources/ADRMPC_outputsVSrefs_inputs.svg)

(a) the $7$ outputs in the ADRMPC simulation with respect to their references values (thin black marked lines); (b) trend of all the $15$ input provided by the ADRMPC during the simulation with respect to their constraints (dotted lines).




![ADRMPC_outputErrorsVSdisturbances.svg](_resources/ADRMPC_outputErrorsVSdisturbances.svg)
Magnitude spectrum of each component of the output error (orange) with respect to the magnitude spectrum of the corresponding component of the additive disturbance (blue) during the simulation.


## How to use the code

* The script `ADRMPCmain10kHz.m` allows to perform the simulation using the ADRMPC framework, exploring the performance for different receding horizon values.
* The script `overactuated_MPC_10kHz.m` allows to perform a comparative simulation, based on an output tracking classical MPC scheme.