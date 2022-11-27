## Fastenating - A Bolted Joint Calculator
Determines factors of safety for fasteners subjected to tensile loading, allows creation of factor of safety plots and joint diagrams. Considers factor of safety against joint separation, tensile yield failure, and fatigue failure. This tool follows the methodology outlined in chapter 15 of the fifth edition of [Norton's Machine Design](https://wps.pearsoned.com/ecs_norton_mechdesign_5/).

### Capability Overview
This calculator considers a single fastener fastening two members where the joint is loaded solely in tension.


<p align="center">
  <img src="https://blog.maxprocorp.com/hs-fs/hubfs/Tension-Load.jpg?width=235&name=Tension-Load.jpg" />
</p>

This calculator is capable of determining factors of safety for varying initial preloads applied to the fastener:
![Preload Diagram](https://github.com/slehmann1/Fastenating/blob/main/SupportingInfo/Images/PreloadDiagram.png?raw=true)
It can also determine factors of safety for the percentage of the fastener proof load that is used to define the initial preload.
![Function of Proof Diagram](https://github.com/slehmann1/Fastenating/blob/main/SupportingInfo/Images/ProofLoadDiagram.png?raw=true)
Finally, the calculator is also capable of plotting standard [bolted joint diagrams](https://www.youtube.com/watch?v=6UowCTm2oT4).
![Bolted Joint Diagram](https://github.com/slehmann1/Fastenating/blob/main/SupportingInfo/Images/JointDiagram.png?raw=true)
#### Determination of Joint Stiffness Factor
The joint stiffness factor is determined through the methodology outlined by [Cornwell](https://journals.sagepub.com/doi/abs/10.1243/09544062JMES1108?journalCode=picb). Instead of the standard conical frustrum approach, this method fits equations to FEA studies of bolts within a variety of plates. In this tool, a key simplifying assumption is made that the members all have equivalent stiffnesses. This may not be valid in all cases, for example where a gasket is introduced into the joint. 

#### Dependencies
Written in python with the following dependencies:  Numpy, Pandas, and MatPlotLib.

#### Limitations
Whilst this tool can be useful for specification of fasteners, all of the assumptions made within it should be understood. An understanding of fastener theory should be grasped before attempting to use this tool; there are multiple factors that can influence fastener failure outside of the parameters this calculator considers. Some of many examples of these factors include temperature, corrosion, or shear, and moment loading. In particular the assumptions made to determine the joint constant and fastener stiffness should be well understood; these assumptions are applicable to select cases and will not perform well in other cases where for example there is a gasket present in the joint. Overall, it would generally be considered poor engineering practice to design bridges based on online and unverified tools.

