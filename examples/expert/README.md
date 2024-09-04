
## Expert

These scenarios demonstrate various examples of complex programs which can be built using pattern discovery or pattern validation functionality offered by Desbordante. These programs aim to provide tangible benefits for end-users by solving real-life problems.

+ [anomaly_detection.py](https://github.com/Desbordante/desbordante-core/tree/main/examples/expert/anomaly_detection.py) — a scenario demonstrating the process of monitoring the quality of newly-arrived data. The idea of this scenario is as follows: first, a user explores a dataset D1 (which is considered a “clean” dataset) by discovering all exact functional dependencies. Next, when the D2 dataset arrives, its exact functional dependencies are discovered and compared with D1’s. There are no missing dependencies, so this portion is considered clean. Then, D3 arrives, the same comparison is performed and the user learns that a single dependency is missing. Afterwards, the user makes various efforts to understand what has happened with data by trying to use approximate functional and, finally, metric functional dependencies. Finally, the user learns that item_id -> item_weight became a metric dependency with range 2. This means that the weight has changed for some items, due to a typo or due to changes in the actual weight. The idea of this approach is sketched in this [paper](https://arxiv.org/abs/2307.14935). There is also a streamlit version of it [here](https://desbordante.streamlit.app/).
+ [dedupe.py](https://github.com/Desbordante/desbordante-core/tree/main/examples/expert/dedupe.py) — a scenario demonstrating how to build a simple entity resolution application using discovered approximate functional dependencies. On each step, the user is presented with a number of candidate table records which they can deduplicate, enrich or delete. The idea of the approach is sketched in these articles ([1](https://itnext.io/building-a-simple-data-cleaning-application-with-desbordante-e4897dcd4c5d#3a03-c6b309624eab), [2](https://arxiv.org/abs/2307.14935)). There is also a streamlit version of it [here](https://desbordante.streamlit.app/).
+ [mine_typos.py](https://github.com/Desbordante/desbordante-core/tree/main/examples/expert/mine_typos.py) — a scenario demonstrating typo discovery inside tables using functional and approximate functional dependencies, which are discovered from a dataset. The user is presented with a number of “almost holding” dependencies, and each of them may indicate the presence of a typo. The script then offers means to easily check them by presenting rows with exceptions (clusters where the exact functional dependency is violated). The idea of this approach is sketched in this [paper](https://arxiv.org/abs/2307.14935). There is also a streamlit version of it [here](https://desbordante.streamlit.app/).