# Phy_Detect_EP
This detection method is proposed in [Expectation Propagation Detection for High-Order High-Dimensional MIMO Systems](https://ieeexplore.ieee.org/document/6841617) by **Javier Céspedes**. It has three process: cavity distribution, Bayesian estimation and iterative update.
> Cespedes, J., Olmos, P. M., Sánchez-Fernández, M., & Perez-Cruz, F. (2014). Expectation propagation detection for high-order high-dimensional MIMO systems. IEEE Transactions on Communications, 62(8), 2840-2849.

Also, the initial version of this code is from `Alva Kosasih`.

## How to install
Currently, we offer three options to install this tool.
* Install through `Matlab Add-Ons`

    Go to [releases](https://github.com/whatshow/Phy_Detect_EP/releases/) to download the file `whatshow_phy_detect_ep.mltbx` in the latest release to install.<br>
    The class `EP` will be automatically added in Matlab. (Don't forget to restart Matlab after installation).
* Install through `pip`
    ```sh
    pip install whatshow-phy-detect-ep
    ```
    * **import this module**
        ```
        from whatshow_phy_detect_ep import EP
        ```
* Install through `git` under another local repositiory
    ```sh
    git submodule add git@github.com:whatshow/Phy_Detect_EP.git Modules/Phy_Detect_EP
    ```
    * **import this module**
        * Matlab
            ```matlab
            addpath("Modules/Phy_Detect_EP");
            ```
        * Python
            ```python
            if '.' not in __name__ :
                from Modules.Phy_Detect_EP.EP import EP
            else:
                from .Modules.Phy_Detect_EP.EP import EP
            ```

## How to use
All EP codes are uniform in matlab and python as a class of `EP`. This class contains the whole process. This section will illustrate the methods of this class following the process of detection.
* EP(the parameter names follow the notations in [Expectation Propagation Detection for High-Order High-Dimensional MIMO Systems](https://ieeexplore.ieee.org/document/6841617))<br>
    `@constellation`: the constellation, a vector<br>
    `@beta`: the percentage for taking values from the previous iteration<br>
    `@epsilon`: the default minimal variance<br>
    `@l`: the maximal iteration<br>
    `@early_stop`: whether stop early<br>
    `@early_stop_min_diff`: the elemental minimal (element-wise) difference in mean and variance<br>
    `@batch_size`**(optional)** : the batch size **(only used in python)**.<br>
    ```matlab
    % matlab
    M = 16;
    constellation = qammod([0: M - 1], M, "UnitAveragePower", true);
    ep = EP(constellation, "beta", 0.1, "epsilon", 1e-13, "early_stop_min_diff", 1e-4);
    ```
    ```python
    # python
    ep = EP(constellation, beta=0.1, epsilon=1e-13, early_stop_min_diff=1e-4); # no batch
    ep = EP(constellation, beta=0.1, epsilon=1e-13, early_stop_min_diff=1e-4, batch_size=10); # using batch
    ```
* detect<br>
    `@y`: the received signal, a vector of [(batch_size), rx_num], [(batch_size), rx_num, 1] or [(batch_size), 1, rx_num] <br>
    `@H`: the channel matrix, a matrix of [(batch_size), rx_num, tx_num]<br>
    `@No`: the noise (linear) power<br>
    `@sym_map`: false by default. If true, the output will be mapped to the constellation<br>
    ```matlab
    # matlab
    syms = ep.detect(y, H, No); # no hard mapping
    syms = ep.detect(y, H, No, sym_map=True); # hard mapping
    ```
    ```python
    # python
    syms = ep.detect(y, H, No); # no hard mapping
    syms = ep.detect(y, H, No, sym_map=True); # hard mapping
    ```

## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require running `init` in the command window first to load directories.
* `Alva`: this is the original code from `Alva Kosasih`.
* `Tests`
    * `Tests/test_case_01`: test the difference between the mature code from and newly proposed `EP` object (no early stop).
    * `Tests/test_case_02`: draw the SER vs SNR between the mature code from and newly proposed `EP` object (early stop).

## Further Reading
### The division of two Gaussian distribution
![latex](https://latex.codecogs.com/gif.latex?\frac{P_1(x)\sim%20N(\mu_1,{\sigma_1}^2)}{P_2(x)\sim%20N(\mu_2,{\sigma_2}^2)})

![latex](https://latex.codecogs.com/gif.latex?=\frac{\frac{1}{\sigma_1\sqrt{2\pi}}e^{-\frac{(x-\mu_1)^2}{2{\sigma_1}^2}}}{\frac{1}{\sigma_2\sqrt{2\pi}}e^{-\frac{(x-\mu_2)^2}{2{\sigma_2}^2}}})

![latex](https://latex.codecogs.com/gif.latex?=\frac{\sigma_2}{\sigma_1}e^{-\frac{1}{2}[\frac{(x-\mu_1)^2}{\sigma_1}%20-\frac{(x-\mu_2)^2}{\sigma_2}]})

Now we set ![latex](https://latex.codecogs.com/gif.latex?X=-\frac{1}{2}[\frac{(x-\mu_1)^2}{\sigma_1}%20-\frac{(x-\mu_2)^2}{\sigma_2}]) and we are going to calculate **`X`** to see the distribution

![latex](https://latex.codecogs.com/gif.latex?X=\frac{1}{\sigma_1^2\sigma_2^2}(\sigma_2^2(x-\mu_1)^2-\sigma_1^2(x-\mu_2)^2))

![latex](https://latex.codecogs.com/gif.latex?=\frac{1}{\sigma_1^2\sigma_2^2}[(\sigma_2^2-\sigma_1^2)x^2-2x(\mu_1\sigma_2^2-\mu_2\sigma_1^2)+(\mu_1^2\sigma_2^2-\mu_2^2\sigma_1^2)])

Therefore, we can conclude that the the new mean and new variance is given as below (ignoring the constant part because this is just the likelihood of two **`Xs`** are same).

![latex](https://latex.codecogs.com/gif.latex?\mu_{new}=\frac{\mu_1\sigma_2^2-\mu_2\sigma1^2}{\sigma_2^2-\sigma_1^2})

![latex](https://latex.codecogs.com/gif.latex?\sigma_{new}^2=(\frac{1}{\sigma_1^2}-\frac{1}{\sigma_2^2})^{-1}=\frac{\sigma_1^2\sigma_2^2}{\sigma_2^2-\sigma_1^2})
### Cavity distribution 
![latex](https://latex.codecogs.com/gif.latex?\hat{p}(y|x)=\frac{p(x|y)}{p(x)}=\frac{N~(x,(H^HH)^{-1}H^Hy,H^H\sigma_w^2H)}{N(x,0,Es)})
### Cavity Marginal
Now we are going to use the **`2 Gaussian Distribution Division`** knowledge to calculate![latex](https://latex.codecogs.com/gif.latex?\hat{p}(y|x)=\frac{p(x|y)}{p(x)}\sim%20N(u_i:t_i,h_i^2))

Now, we follow the Math symbols in the paper, suppose  ![latex](https://latex.codecogs.com/gif.latex?p(x|y)\sim%20N(\mu_i,\sigma_i^2)) ( ![latex](https://latex.codecogs.com/gif.latex?\sigma_i^2) is the diagonal of the matrix ![latex](https://latex.codecogs.com/gif.latex?\Sigma) ) and  ![latex](https://latex.codecogs.com/gif.latex?p(x)\sim%20N(\gamma\Lambda^{-1},\Lambda^{-1})) 

Here, we need pay attention that the mean of **`p(x)`** is ![latex](https://latex.codecogs.com/gif.latex?\gamma\Lambda^{-1}) instead of ![latex](https://latex.codecogs.com/gif.latex?\gamma). The proof is give below where we follow the exponential family.

![latex](https://latex.codecogs.com/gif.latex?e^{\gamma_i%20u_i-\frac{1}{2}\Lambda_i%20u_i^2}=e^{-\frac{u_i^2-2\gamma\Lambda^{-1}%20u_i}{2\Lambda^{-1}}})

Here, we can see that the mean is ![latex](https://latex.codecogs.com/gif.latex?\gamma\Lambda^{-1})

Now,  we are going to find out  ![latex](https://latex.codecogs.com/gif.latex?N(u_i:t_i,h_i^2)) in *`Lth`* iteration

![latex](https://latex.codecogs.com/gif.latex?h_i^2=\frac{\sigma_1^2\sigma_2^2}{\sigma_2^2-\sigma_1^2}=(\frac{1}{\sigma_i^2}-\frac{1}{\Lambda^{-1}})^{-1}=(\frac{1}{\sigma_i^2}-\Lambda)^{-1}=(\frac{1}{\sigma_i^2}-\frac{\sigma_i^2\Lambda}{\sigma_i^2})^{-1}=\frac{\sigma_i^2}{1-\sigma_i^2\Lambda})

![latex](https://latex.codecogs.com/gif.latex?t_i=\frac{\mu_1\sigma_2^2-\mu_2\sigma_1^2}{\sigma_2^2-\sigma_1^2}=\frac{\sigma_1^2\sigma_2^2}{\sigma_2^2-\sigma_1^2}\frac{\mu_1\sigma_2^2-\mu_2\sigma_1^2}{\sigma_2^2\sigma_1^2}=h_i^2\frac{\mu_1\sigma_2^2-\mu_2\sigma_1^2}{\sigma_1^2\sigma_2^2})

![latex](https://latex.codecogs.com/gif.latex?=h_i^2(\frac{\mu_1}{\sigma_1^2}-\frac{\mu_2}{\sigma_2^2})=h_i^2(\frac{\mu_i}{\sigma_i^2}-\frac{\gamma\Lambda^{-1}}{\Lambda^{-1}})=h_i^2(\frac{\mu_i}{\sigma_i^2}-\gamma))
### The Pair Update ![latex](https://latex.codecogs.com/gif.latex?(\gamma_i^{l+1},\Lambda_i^{l+1}))
From the step computing the mean ![latex](https://latex.codecogs.com/gif.latex?\hat{\mu}_i) and the variance ![latex](https://latex.codecogs.com/gif.latex?\hat{\sigma}_i^2) of prediction ![latex](https://latex.codecogs.com/gif.latex?\hat{p}(x|y)), we can use those 2 results to calculate ![latex](https://latex.codecogs.com/gif.latex?p(x)^{l+1}); that is, for *`(L+1)th`* iteration,

![latex](https://latex.codecogs.com/gif.latex?p(x)^{(l+1)}\sim%20N(\frac{\gamma_i^{(l+1)}}{\Lambda_i^{(l+1)}},\frac{1}{\Lambda_i^{(l+1)}})=\frac{\hat{p}(x|y)^{(l)}\sim%20N(\mu_{pi},\sigma_{pi}^2)^{(l)}}{p(y|x)^{(l)}\sim%20N(t_i,h_i^2)^{(l)}})

Using  **`2 Gaussian Distribution Division`** knowledge, for the variance, 

![latex](https://latex.codecogs.com/gif.latex?\frac{1}{\Lambda_i^{(l+1)}}=(\frac{1}{\sigma_1^2}-\frac{1}{\sigma_2^2})^{-1})

![latex](https://latex.codecogs.com/gif.latex?\Lambda_i^{(l+1)}=\frac{1}{\sigma_1^2}-\frac{1}{\sigma_2^2}=\frac{\sigma_2^2-\sigma_1^2}{\sigma_1^2\sigma_2^2}=\frac{1}{\sigma_{pi}^{2(l)}}-\frac{1}{h_i^{2(l)}})
For the mean,

![latex](https://latex.codecogs.com/gif.latex?\frac{\gamma_i^{(l+1)}}{\Lambda_i^{(l+1)}}=\frac{\mu_1\sigma_2^2-\mu_2\sigma1^2}{\sigma_2^2-\sigma_1^2})

![latex](https://latex.codecogs.com/gif.latex?\gamma_i^{(l+1)}=\Lambda_i^{(l+1)}\frac{\mu_1\sigma_2^2-\mu_2\sigma1^2}{\sigma_2^2-\sigma_1^2}=\frac{\sigma_2^2-\sigma_1^2}{\sigma_1^2\sigma_2^2}\frac{\mu_1\sigma_2^2-\mu_2\sigma1^2}{\sigma_2^2-\sigma_1^2}=\frac{\mu_1\sigma_2^2-\mu_2\sigma1^2}{\sigma_1^2\sigma_2^2})

![latex](https://latex.codecogs.com/gif.latex?=\frac{\mu_1}{\sigma_1^2}-\frac{\mu_2}{\sigma_2^2}=\frac{\mu_{pi}^{(l)}}{\sigma_{pi}^{2(l)}}-\frac{t_i^{(l)}}{h_i^{2(l)}})