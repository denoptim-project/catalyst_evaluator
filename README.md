# Automated Evaluation of Ru-Based Catalyst for Olefin metathesis
This is a tool for automatic evaluation of candidate olefin metathesis catalyst with general formula (L)Ru(Y)(X)=CH<sub>2</sub>, where X and Y are covalent ligands and L is a dative ligand.

TODO: add figure with workflow overview

## How To Get Started
1. Clone/copy this repository on your local client.
    ```
    git clone --recurse-submodules <git_or_https_URL>
    ```

2. The file [environment.yml](environment.yml) can be used to create a suitable running environment. This is best done with [mamba](https://mamba.readthedocs.io/en/latest/)
    ```
    mamba env create --file environment.yml
    ```
    but can also be done with [conda](https://docs.conda.io/en/latest/index.html):
    ```
    conda env create --file environment.yml
    ```
    After creation, you should activate the environment with
    ```
    conda activate RuCatEvaluator
    ```

3. In addition to the pre-configured environment, <a href="https://dasher.wustl.edu/tinker/">Tinker</a> and <a href="https://www.wavefun.com/spartan">Spartan</a> need to be installed according to the corresponding vendor as neither of them can be installed via Conda. Please refer to the above links for license terms and installation instructions. After installation, make the executables reachable, by adjusting the corresponding values in <a href="evaluate_candidate.sh">evaluate_candidate.sh</a>:
    ```
    export TINKERBIN="_your_path_to_Tinkers_bin_folder/bin"
    export SPARTANEXE="_your_path_to_Spartan_executable/spartan20"
    ```

4. We now configure some remote computer to run Gaussian DFT and, possibly, xTB calculations. These remote computers are typically HPCs. The interface with the remote computers is managed by the tool [RemoteWorkersBridge](https://github.com/denoptim-project/RemoteWorkersBridge) (or git submodule RemoteWorkersBridge under the [tools](tools) folder). In general, we assume you have a way to send jobs to a queuing system on such HPC workers or start such jobs in whichever way according to what is suitable for your specific remote worker. In our case this task is performed by the command `submit_job_acc`, which we assume you'll make available in your HPC workers (TODO: make source available at [tools/submit_job_acc/submit_job_acc.sh](tools/submit_job_acc/submit_job_acc.sh)). The following steps allow to configure the bridge to a remote worker (in alternative, have a look at the [Test Run Without Remote Workers](#test-run-without-remote-workers)).

5. Create a pair of ssh keys
    ```
    ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa_RuCatEvaluator
    ```

6. Add the key to an agent:
    ```
    eval `ssh-agent -s`
    ssh-add ~/.ssh/id_rsa_RuCatEvaluator
    ```

7. Identify the IP of the remote worker, declare your user name on that machine, and specify a pathname that can be used to place job files on the remote worker. Edit the following commands and run them on your local machine:
   ```
   export MYRWIP=<your_worker_IP>
   export MYWUSER=<your_username>
   export MYWRKDIR=<path_to_remote_work_space>
   ```

7. Enable the key on the remote host:
    ```
    ssh-copy-id -i ~/.ssh/id_rsa_RuCatEvaluator $MYWUSER@$MYRWIP
    ```
    where `your_username` and `your_worker_IP` should be replaced with your specific user-name and IP address. Presently we support only IPv4. You can get the proper IP by running `echo $(curl -s -4 ifconfig.me/ip)` on the remote.

8. Configure the bridge. Run the following command in your local machine. Set `MYWRKDIR` to an absolute pathname to the folder located on the remote worker: this is the location from where every job on the remote worker will start.


    ```
    ssh -i ~/.ssh/id_rsa_RuCatEvaluator $MYRWIP  "mkdir -p $MYWRKDIR"
    scp -r -i ~/.ssh/to_RuCatEvaluator data/basisset MYRWIP:"$MYWRKDIR"
    ```
    and
    ```
    echo "[WORKER1]" > tools/RemoteWorkersBridge/configuration
    echo "remoteIP=your_worker_IP" >> tools/RemoteWorkersBridge/configuration
    echo "wdirOnRemote=$MYWRKDIR" >> tools/RemoteWorkersBridge/configuration
    echo "userOnRemote=$MYWUSER" >> tools/RemoteWorkersBridge/configuration
    echo "identityFile=$HOME/.ssh/id_rsa_RuCatEvaluator" >> tools/RemoteWorkersBridge/configuration
    echo "workKind=xtb,dft" >> tools/RemoteWorkersBridge/configuration
    ```

9. Copy the folder `tools/RemoteWorkersBridge` into a private location on the remote worker, typically under your HOME. We'll refer to the absolute path of the newly created folder as `RemoteWorkersBridge_on_HPCWorker`, which you have to edit in this command:
    ```
    scp -r -i ~/.ssh/to_RuCatEvaluator tools/RemoteWorkersBridge MYRWIP:<RemoteWorkersBridge_on_HPCWorker>
    ```

10. Log in to `your_worker_IP` and edit the `~/.ssh/authorized_keys` file. The last line of this file should contain the ssh key entry you have just added with the `ssh-copy-id` command above. We are now going to edit this line to prevent any misuse of this automated login channel. This is done by limiting the usage to this key enabling only a privately own command filter. To this end, edit the line pertaining the ssh key we just authorized (i.e., the last line of `~/.ssh/authorized_keys`), and add in front of any text of that line the following string (NB: there is a space at the end!):
    ```
    from="your_IP",command="RemoteWorkersBridge_on_HPCWorker/commandFilter.sh"
    ```
    where `your_IP` is the IP address of your local client (the machine what will use this connection to submit jobs to the worker. Presently we support only IPv4. You can get the proper IP by running `echo $(curl -s -4 ifconfig.me/ip)` from within your local client) and `RemoteWorkersBridge_on_HPCWorker` is the absolute path defined in the previous step.

11. Test the interface with the remote worker from your local client:
    ```
    cd tools/RemoteWorkersBridge/submit_tool/test/
    ./runTest.sh
    ```

12. Depending on the needs, you may want to configure the scripts to check for completion of remote jobs with high frequency (useful only when testing if everything is working fine) and if the xTB jobs should be run locally or in the remote worker. These settings can be controlled in the initial part of the `evaluate_candidate.sh` script.

13. You should be ready to go now!


## Test Run
A pre-configured test run is available. The time required to run it depends the loading on the remote workers, but it should be not shorted than a couple of hours. So, be patient and make sure you can keep the session alive (you may want to look into tools that allow to detach the terminal session, e.g., [tmux](https://github.com/tmux/tmux/wiki):
```
./test/run_test.sh
```    
This with run for a few instants and then return `Test PASSED!` if the entire workflow could be successfully executed and the resulting fitness is sufficiently close to the expected value.


## Test Run Without Remote Workers
This section is for those who, driven by testing and/or development needs, may want to run the evaluator without configuring a proper remote worker and without really running any DFT job. To this end, we can use the local client, i.e., `localhost` as a fake remote worker. We can build an ssh bridge to such fake remote machine and ask it to deliver results for a pre-configured test run.

Here are the steps to follow:
1. Ensure we can login on `localhost`. We need to enable remote login __and__ to configure an ssh key.
    - The enabling remote login may imply adjusting the firewall settings or explicitly enabling this possibility. On many machine that you routinely login to via ssh, this is most likely not a problem. However, this is typically needed on standard laptops. How to do this, it depends on your OS. For example, MacOS you need to do _System Preferences_ -> _Sharing_ -> _Remote Login_ to select which user can login via ssh.
    - The creation of an ssh key can be done by `ssh-keygen -t rsa -b 4096 -f ~/.ssh/to_localhost` specifying an empty pass-phrase, followed by `cat ~/.ssh/to_localhost.pub >> ~/.ssh/authorized_keys`. At this point, the following should allow you to login: `ssh -i ~/.ssh/to_localhost localhost`. Logout again and proceed.

2. From now on, we'll proceed as in the normal setup but we treat the localhost as a fake remote worker.

3. Define a workspace to be created on the fake remote worker. Edit the `<path_to_fake_remote_work_space>` part and run the following from within the folder where this README file is located:
    ```
    export MYWRKDIR=<path_to_fake_remote_work_space>
    ssh -i ~/.ssh/to_localhost localhost  "mkdir -p $MYWRKDIR"
    scp -r -i ~/.ssh/to_localhost data/basisset localhost:"$MYWRKDIR"
    ```

4. Configure the ssh bridge to the fake remote worker:
    ```
    echo "[WORKER1]" > tools/RemoteWorkersBridge/configuration
    echo "remoteIP=localhost" >> tools/RemoteWorkersBridge/configuration
    echo "wdirOnRemote=$MYWRKDIR" >> tools/RemoteWorkersBridge/configuration
    echo "userOnRemote=$USER" >> tools/RemoteWorkersBridge/configuration
    echo "identityFile=$HOME/.ssh/to_localhost" >> tools/RemoteWorkersBridge/configuration
    echo "workKind=xtb,dft" >> tools/RemoteWorkersBridge/configuration
    ```

5. Place a copy of the customised `RemoteWorkersBridge` folder into the fake remote worker:
    ```
    scp -r -i ~/.ssh/to_localhost tools/RemoteWorkersBridge localhost:"$HOME/RemoteWorkersBridge_to_localhost"
    ```

6. Remove the last line of `~/.ssh/authorized_keys` (__NB: we assume no other change has occurred to this file since point 1 in this workflow__).

7. Finally, enforce the use of a command filter to ssh commands coming via the ssh bridge:
    ```
    echo "command=\"$HOME/RemoteWorkersBridge_to_localhost/commandFilter.sh\" $(cat ~/.ssh/to_localhost.pub)" >> ~/.ssh/authorized_keys
    ```

8. Optionally, we may want to verify the functionality of the newly-created ssh bridge:
    ```
    cd tools/RemoteWorkersBridge/submit_tool/test
    ./runTest.sh
    cd ../../../..
    ```
      This will print a comforting message if the test has been passed.

9. To make the localhost emulate an HPC worker we need to define the commands that run or submit computational chemistry jobs on the HPC resource. Such commands are left to the user, which will have to create them according to the specifics of the HPC resource. Such commands are `submit_job_acc-2`, which performs AutoCompChem-controlled multistep workflows by sending it to a job scheduler, and `run_xtb`, which runs the xTB software without sending it to a job scheduler. Notably, all the `run_*` commands are meant to start a calculation on the cpus that are available to the process that calls them. Overall, the strategy is to add the corresponding executables to the PATH to make them reachable. This is here done by altering the ~/.bashrc file (assuming BASH terminal).  __NB: this step makes this workflow incompatible with situation where you want to have alternative versions of such commands. So, you may want to recover the original ~/.bashrc file once you are done with this test.__ Note that depending on your OS, your ~/.bashrc may exclude non-interactive shells. In such case, this step will fail unless you make the ~/.bashrc be effective also on non-interactive shells. Also, note we assume we are in the root folder of the repository, i.e., right where this README file is located.
    ```
    echo "# Added to use localhost to emulate HPC worker:" >> ~/.bashrc
    echo "export PATH=\"$(pwd)/tools/hpc_emulator:$(pwd)/tools/xtb_runner:\$PATH\"" >> ~/.bashrc
    ```

10. Now we start the test run as if we were submitting to a remote worker:
    ```
    ./test/run_test.sh
    ```
      This with run for a few instants and then return `Test PASSED!` if the entire workflow could be successfully executed.


## Evaluation of Catalysts
__NB.__ The time required to evaluate a single catalyst depends on the loading on the remote workers, but is typically of few hours. To make sure you can keep the terminal session alive you may want to look into tools that allow to detach the terminal session, e.g., [tmux](https://github.com/tmux/tmux/wiki).

To define a catalyst with general formula (L)Ru(Y)(X)=CH<sub>2</sub>, where X and Y are covalent ligands and L is a dative ligand, we use [Denoptim's graph representation](https://denoptim-project.github.io/DENOPTIM/) (there is an introductory [tutorial here](https://denoptim-project.github.io/tutorials/tutorial_1.1.html)). Such representation can be prepared in two different ways:

- manually assembling fragments.
    1. Lauch Denoptim's GUI:
    ```
    cd data
    denoptim
    ```
    2. Choose `File` -> `New Graphs` and click on `Load BBSpace`.
    3. Select `Use parameters from existing file` and navigate to load file `data/FS-C5.0.par`.
    4. Now click on `Add Graphs`, select `Build` and start from a `scaffold`. In the dialog click on `Select current vertex` to initiate the graph.
    5. SHIFT-click on the yellow point by the "AP1" and "AP2" labels to select the attachment points meant for X-type ligands, click `Add vertex from BBSpace`, then choose `Compatible vertices`. Navigate to chose the X ligands and confirm your selection.
    6. Select the yellow point that represents "AP3", click `Add vertex from BBSpace`, and then choose `Compatible vertices`. Proceed repeating this step until the L-type ligand of your choice is complete.
    7. Finally, `Export graph` in `SDF` format as `my_catalyst.sdf`. This file can be fed to the catalyst evaluator script (see below).

- chopping a molecular model.


Once you have a Denoptim graph file (e.g., `my_catalyst.sdf`) you can start the evaluation with the following command:
```
./evaluate_candidate.sh -i my_catalyst.sdf
```


# Acknowledgements
The Research Council of Norway and University of Bergen are thanked for various kinds of funding.
