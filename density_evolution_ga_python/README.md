# Background


An optimal LDPC code is typically defined by a waterfall region close to the Shannon limit of the channel. In principle, when finding an optimal code that performs well in the waterfall region, it would be possible to simulate a large number of codes and pick the code that performs the best. However, due to the large search space and the time it takes to simulate, this would take an immensely long time. Therefore, it is necessary to find a method in which the decoding performance at the waterfall region can be characterized through probabilistic measures of groups of codes. By looking at code ensembles, the search space becomes smaller and simulation times become much faster, which makes the optimization computationally feasible. To further motivate the use of code ensembles, it can be shown that decoding failures that characterize the type of failuresthat occur at the waterfall region are sharply concentrated around the ensemble average. This makes code ensembles representative of the decoding performance at the waterfall region of individual codes within.


# Program description
The main script is [density_evolution.py](./density_evolution.py). The parameters that are meant to be adjusted by the user (along with a parameter description) can be found in [config.yml](./config.yml). Furthermore, [de_utils.py](./de_utils.py) contains code that implements density evolution along with some other functions that are used for optimization.

There are two types of algorithms implemented to optimize degree distributions based on density evolution. The first ([ga_continuous.py](./ga_continuous.py)) uses the concept of differential evolution to find optimal solutions constrained by maximum allowed degree distributions. Furthermore, the degree distributions are considered continuous in this aspect. This algorithm is not error free and could be improved significantly.  

The other algorithm, and the one that has been used most in this project ([ga_discrete.py](ga_discrete.py)), genetically optimizes degree distributions based on the constraints set by the size of the protograph used in [MM_QC_PEGA](../mm_qc_pega).

Lastly, since the algorithms are rather complex, the multiprocessing library has been used to run the optimizations in parallel. An illustration of how multiprocessing has been used can be seen in [multiprocessing_illustration.py](multiprocessing_illustration.py).

## Dealing with asymmetry of the channel
It is important to mention that the concept used in this project to deal with asymmetry has been the one described on page 312 in [Modern Coding Theory](https://books.google.se/books?hl=en&lr=&id=ZJrZPObOe60C&oi=fnd&pg=PR13&dq=modern+coding+theory+urbanke&ots=WohuOu5lqr&sig=OaewAYIeWVBUZQYBnexWfgaHo3k&redir_esc=y#v=onepage&q&f=false), i.e. the channel has been symmetrized and then evaluated using symmetrical density evolution. This yielded some unexpected results, and it would be interesting to see that the results would be the same for asymmetric density evolution as described in [Wang's article](https://ieeexplore.ieee.org/abstract/document/1542413).

# Running scripts on a remote server
## Install Anaconda on remote server
To be able to run the script on a remote server you will need the correct version of python. In our versions, we have used python 3.8.5. Installing directly on the OS is not adviced, so install using some virtual environment, for example Anaconda. Here, we go trhough how to install with Anaconda.

Install using wget:
 - Go to [Anaconda archive](https://repo.anaconda.com/archive/) and pick the latest version of Anaconda compatible with linux (Linux-x86_64.sh)
 - Log into your remote server using ssh.
 - In your root folder, type:
  ```
  $ wget https://repo.anaconda.com/archive/[Anaconda_version.sh]
  ```

 - Replace Anaconda_version.sh with the one you picked from the archive.
 - Next install Anaconda (executing in your root folder):
  ```
  $ bash [Anaconda_version.sh]
  ```

 - Accept the license and confirm installation location.
 - Next, you will need to set you conda environment variable. Do this by typing
  ```
  $ echo 'export PATH=~/anaconda3/bin:$PATH' >> ~/.bashrc
  $ source ~/.bashrc
  ```

You should now be all set to start using Anaconda. Check documentation for how to run different commands. 

Essentiallly, you want to create an environment for example

```
$ conda create --name ldpc_master_thesis python=3.8.5
$ source activate ldpc_master_thesis
```

Check that you are in the correct environment and start installing using conda install. You should now have the correct version of python when running scripts.



## Run jobs using screen

To run script remotely, use "screen" on linux.
To start a new screen type:

```
$ screen -S [screen_name]
```

You can pick your own "screen_name". 

Then, start your job inside the screen, and if you want to keep it running even after you have logged out from ssh press Ctrl+a+d.

To go into your screen again, type

```
$ screen -R [screen_name]
``` 

If you have forgotten your screen name:


```
$ screen -ls
```

# Sources

Design of capacity-approaching irregular low-density parity-check codes: [https://ieeexplore.ieee.org/document/910578](https://ieeexplore.ieee.org/document/910578)

Density Evolution for Asymmetric Memoryless Channels: [https://arxiv.org/pdf/cs/0509014.pdf](https://arxiv.org/pdf/cs/0509014.pdf)
