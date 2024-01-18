# Automated Evaluation of Ru-Based Catalyst for Olefin metathesis
This is a tool for automatic evaluation of candidate olefin metathesis catalyst with general formula (L)Ru(Y)(X)=CH<sub>2</sub>, where X and Y are covalent ligands and L is a dative ligand.

TODO: add figure with workflow overview

## How To Get Started
1. Clone/copy this repository on your local client.
    ```
    git clone --recurse-submodules git@github.com:denoptim-project/catalyst_evaluator.git
    ```
    and move into it
    ```
    cd catalyst_evaluator
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

3. In addition to the pre-configured environment, <a href="https://dasher.wustl.edu/tinker/">Tinker</a> and <a href="https://www.wavefun.com/spartan">Spartan</a> need to be installed according to the corresponding vendor as neither of them can be installed via Conda. Please refer to the above links for license terms and installation instructions. After installation, make the executables reachable, by adjusting the corresponding values in <a href="evaluate_catalyst.sh">evaluate_catalyst.sh</a>:
    ```
    export TINKERBIN="_your_path_to_Tinkers_bin_folder/bin"
    export SPARTANEXE="_your_path_to_Spartan_executable/spartan20"
    ```

4. We now configure a remote computer (later referred to as a "remote worker") to run Gaussian DFT and, possibly, xTB calculations, which are typically executed in HPC service provider. The interface with one or more remote computers is managed by the tool [RemoteWorkersBridge](https://github.com/denoptim-project/RemoteWorkersBridge) (or git submodule RemoteWorkersBridge under the [tools](tools) folder). In general, we assume you have a way to send jobs to a queuing system on such HPC workers or start such jobs in whichever way according to what is suitable for your specific remote worker. In our case this task is performed by the command `submit_job_acc`, which we assume you'll make available in your HPC workers (TODO: make source available at [tools/submit_job_acc/submit_job_acc.sh](tools/submit_job_acc/submit_job_acc.sh)). The following steps allow to configure the bridge to a remote worker (in alternative, have a look at the [Test Run Without Remote Workers](#test-run-without-remote-workers) ).

5. Create two pairs of ssh keys (with non-empty pass-phrase):
    ```
    ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa_RuCatEvaluator
    ```
    and
    ```
    ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa_RuCatEvaluator_setup
    ```

6. Add both keys to an agent:
    ```
    eval `ssh-agent -s`
    ssh-add ~/.ssh/id_rsa_RuCatEvaluator
    ```
    and
    ```
    ssh-add ~/.ssh/id_rsa_RuCatEvaluator_setup
    ```

7. Configure a remote worker by defining the following variables. Replace each `<to_be_replaced>` with the appropriate string and run the resulting command in your local machine.
    - The IP of the remote worker. Presently we support only IPv4. You can get the proper IP by running `echo $(curl -s -4 ifconfig.me/ip)` on the remote worker.
    ```
    export MYWORKERIP=<to_be_replaced>
    ```
    - The user name to use to connect to the remote worker. This may or may not be different from your user name of the local machine.
    ```
    export MYWORKERUSER=<to_be_replaced>
    ```
    - The absolute pathname of a non-existing folder on the remote worker. this is where we'll send all input files for the jobs to be run by the worker.
    ```
    export MYWORKERDIR=<to_be_replaced>
    ```
    - The absolute pathname of a non-existing folder for hosting the remote end of the RemoteWorkersBridge tool. This is typically under the HOME or a location that is not shared with other users.
    ```
    export MYWORKERBRIDGEHEAD=<to_be_replaced>
    ```

8. Authorize the setup key on the remote worker:
    ```
    ssh-copy-id -i ~/.ssh/id_rsa_RuCatEvaluator_setup $MYWORKERUSER@$MYWORKERIP
    ```

9. Complete the configuration by running the following from the base folder of this repository in your local machine.
    ```
    echo "[WORKER1]" > tools/RemoteWorkersBridge/configuration
    echo "remoteIP=$MYWORKERIP" >> tools/RemoteWorkersBridge/configuration
    echo "wdirOnRemote=$MYWORKERDIR" >> tools/RemoteWorkersBridge/configuration
    echo "userOnRemote=$MYWORKERUSER" >> tools/RemoteWorkersBridge/configuration
    echo "identityFile=$HOME/.ssh/id_rsa_RuCatEvaluator" >> tools/RemoteWorkersBridge/configuration
    echo "workKind=xtb,dft" >> tools/RemoteWorkersBridge/configuration
    ```

10. To prepare the remote worker to accept incoming requests to run selected commands (see [the command filter](tools/RemoteWorkersBridge/commandFilter.sh) ), run the following five commands from your local machine:
    ```
    ssh -i ~/.ssh/id_rsa_RuCatEvaluator_setup $MYWORKERUSER@$MYWORKERIP  "mkdir -p $MYWORKERDIR"
    ```
    ```
    scp -r -i ~/.ssh/id_rsa_RuCatEvaluator_setup data/basisset $MYWORKERUSER@$MYWORKERIP:"$MYWORKERDIR"
    ```
    ```
    scp -r -i ~/.ssh/id_rsa_RuCatEvaluator_setup tools/RemoteWorkersBridge $MYWORKERUSER@$MYWORKERIP:"$MYWORKERBRIDGEHEAD"
    ```
    ```
    ssh -i ~/.ssh/id_rsa_RuCatEvaluator_setup $MYWORKERUSER@$MYWORKERIP  "chmod go-rwx $MYWORKERBRIDGEHEAD"
    ```
    ```
    echo "from=\"$(curl -s -4 ifconfig.me/ip)\",command=\"$MYWORKERBRIDGEHEAD/commandFilter.sh\" $(cat ~/.ssh/id_rsa_RuCatEvaluator.pub)" | ssh -i ~/.ssh/id_rsa_RuCatEvaluator_setup $MYWORKERUSER@$MYWORKERIP "cat >> ~/.ssh/authorized_keys"
    ```

11. Test the interface with the remote worker from your local client:
    ```
    ./tools/RemoteWorkersBridge/submit_tool/test/runTest.sh
    ```
    This should terminate with a comforting "Test PASSED!" message.

12. Depending on your machine and needs, you may want to set additional parameters in the `evaluate_catalyst.sh` script.

13. You should be ready to go now!


## Test Run
A pre-configured test run is available. The time required to run it depends the loading on the remote workers, but it should be not shorted than a couple of hours. So, be patient and make sure you can keep the session alive (you may want to look into tools that allow to detach the terminal session, e.g., [tmux](https://github.com/tmux/tmux/wiki):
```
./test/run_test.sh
```    
This eventually return `Test PASSED!` if the entire workflow could be successfully executed and the resulting fitness is sufficiently close to the expected value.


## Test Run Without Remote Workers
This section is for those who, driven by testing and/or development needs, may want to run the evaluator without configuring a proper remote worker and without really running any DFT jobs. To this end, we can use the local client, i.e., `localhost` as a fake remote worker. We can build an ssh bridge to such fake remote machine and ask it to deliver results for a pre-configured test run.

Here are the steps to follow:
1. Ensure you can ssh to `localhost`. We need to enable remote login __and__ to configure an ssh key.
    - The enabling remote login may imply adjusting the firewall settings or explicitly enabling this possibility in your machine's configuration. On many machine that you routinely login to via ssh, this is most likely not a problem. However, this is typically needed on standard laptops. How to do this, it depends on your OS. For example, MacOS you need to do _System Preferences_ -> _Sharing_ -> _Remote Login_ to select which user can login via ssh.
    - The creation of an ssh key can be done by the following (NB: specifying an empty pass-phrase)
      ```
      ssh-keygen -t rsa -b 4096 -f ~/.ssh/to_localhost
      ```

2. Now, the configuration of the remote bridge is straightforward: run the following from the base folder of this repository in your local machine to set it up.
    ```
    echo "[WORKER1]" > tools/RemoteWorkersBridge/configuration
    echo "remoteIP=localhost" >> tools/RemoteWorkersBridge/configuration
    echo "wdirOnRemote=$(pwd)/test/localhost_workspace" >> tools/RemoteWorkersBridge/configuration
    echo "userOnRemote=$USER" >> tools/RemoteWorkersBridge/configuration
    echo "identityFile=$HOME/.ssh/to_localhost" >> tools/RemoteWorkersBridge/configuration
    echo "workKind=xtb,dft" >> tools/RemoteWorkersBridge/configuration
    ```

3. Prepare the folder structure
    ```
    mkdir -p test/localhost_workspace
    cp -r data/basisset test/localhost_workspace
    ```
    and force the interpretation of commands from sent via the remote workers RemoteWorkersBridge
    ```
    echo "command=\"$(pwd)/tools/RemoteWorkersBridge/commandFilter.sh\" $(cat ~/.ssh/to_localhost.pub)" >> ~/.ssh/authorized_keys
    ```

4. Test the interface with the localhost via the remote workers bridge:
    ```
    ./tools/RemoteWorkersBridge/submit_tool/test/runTest.sh
    ```
    This will print a comforting message if the test has been passed.

5. To make the localhost emulate an HPC worker we need to define the `submit_job_acc` command that we expect to find in any remote worker. Therefore, we add the corresponding executables to the PATH to make them reachable. This is here done by altering the `~/.bashrc` file (assuming BASH terminal).  __NB: this step makes this workflow incompatible with situation where you want to have alternative versions of such commands. So, you may want to recover the original `~/.bashrc` file once you are done with this test.__ Because of the usual practice of terminating `~/.bashrc` prematurely for non-interactive shells (see [here](https://unix.stackexchange.com/questions/257571/why-does-bashrc-check-whether-the-current-shell-is-interactive) ), you __must__ make sure you add the following in `~/.bashrc` before the line where it checks whether the shell is interactive or not. Also, note you __must__ replace `<replace_with_path>` with the appropriate pathname, i.e., the absolute pathname of the root folder of this repository, i.e., right where this README file is located.
    ```
    export PATH="<replace_with_path>/tools/hpc_emulator:<replace_with_path>/tools/xtb_runner:$PATH"
    ```

6. Now we start the test run as if we were submitting to a remote worker:
    ```
    ./test/run_test.sh --sendXtbToRemote --highFreq
    ```
    This with run for a few instants and then return `Test PASSED!` if the entire workflow could be successfully executed. However, this configuration without a true remote worker is only available for the test run: it cannot be used to evaluate any other catalyst.

If you need to change between running tests with actual HPC workers or localhost, and you have configured both of them, you only need to replace the `tools/RemoteWorkersBridge/configuration` file.


## Evaluation of Catalysts
__NB.__ The time required to evaluate a single catalyst depends on the loading on the remote workers, but is typically of few hours. To make sure you can keep the terminal session alive you may want to look into tools that allow to detach the terminal session, e.g., [tmux](https://github.com/tmux/tmux/wiki).

To define a catalyst with general formula (L)Ru(Y)(X)=CH<sub>2</sub>, where X and Y are covalent ligands and L is a dative ligand, we use [Denoptim's graph representation](https://denoptim-project.github.io/DENOPTIM/) (there is an introductory [tutorial here](https://denoptim-project.github.io/tutorials/tutorial_1.1.html)). Such representation can be prepared in two different ways:

- manually assembling fragments.
    1. Lauch Denoptim's GUI:
    ```
    cd data
    denoptim
    ```
    2. Choose `File` -> `New` -> `New Graphs` and click on `Load BBSpace`.
    3. Select `Use parameters from existing file` and navigate to load file `data/FS-C5.0.par`.
    4. Now click on `Add Graphs`, select `Build` and start from a `scaffold`. In the dialog click on `Select current vertex` to initiate the graph.
    5. SHIFT-click on the yellow point by the "AP1" and "AP2" labels to select the attachment points meant for X-type ligands, click `Add vertex from BBSpace`, then choose `Compatible vertices`. Navigate to chose the X ligands and confirm your selection.
    6. Select the yellow point that represents "AP3", click `Add vertex from BBSpace`, and then choose `Compatible vertices`. Proceed repeating this step until the L-type ligand of your choice is complete.
    7. Finally, `Export graph` in `SDF` format as `mycatalyst.sdf` (or any other filename, but avoid using character "_" as it is later used as separator). This file can be fed to the catalyst evaluator script (see below).

- chopping an existing molecular model of a catalyst with general formula (L)Ru(Y)(X)=CH<sub>2</sub> or the 2-isopropoxybenzylidene precursor (e.g., the classical Hoveyda-Grubbs second generation catalyst).
    1. Lauch Denoptim's GUI:
    ```
    cd data
    denoptim
    ```
    2. Choose `File` -> `New` -> `New Graphs` and click on `Add Graphs` and then `Convert` to select which molecular structure to load. __NB:__ such model is assumed to adhere to the tradition of reporting any Ru-L and Ru-X bond as single bonds and the Ru-alkylidene bond as a double one.
    3. In the resulting dialog window, click on `Import Rules` and navigate to load the rules from the file given in `data/cutting_rules_for_mol-to-graph`.
    4. Specify the policy for identifying the scaffold vertex to `ELEMENT` and specify elemental symbol `Ru`.
    5. By clicking on `Start Fragmentation` to convert the molecule into a graph, i.e., a collection of building blocks.
    6. Make sure that the resulting graph contains three ligands:
        - a scaffold containing Ru,
        - two fragments connected to the scaffold via a `RuX:0-RuX:1` AP-AP connection, and
        - one fragment connected to the scaffold via a `RuL:0-RuL:1` AP-AP connection.
    7. Finally, `Export graph` in `SDF` format as `mycatalyst.sdf` (or any other filename, but avoid using character "_" as it is later used as separator). This file can be fed to the catalyst evaluator script (see below).


Once you have a Denoptim graph file, i.e., the `mycatalyst.sdf` mentioned above, you can start the evaluation with the following command:
```
./evaluate_catalyst.sh -i mycatalyst.sdf
```
To see how the master job proceeds, check the log file `mycatalyst/mycatalyst_FProvider.log`. For a successfully terminated evaluation, the final results are summarised in file `mycatalyst/mycatalyst_out.sdf` while the details are collected in archive `mycatalyst/mycatalyst.tar.gz`.


# Acknowledgements
The Research Council of Norway and University of Bergen are thanked for various kinds of funding.
