# maxflow
Efficient C implementation of the Ford-Fulkerson maxflow algorithm


# Info:
Feel free to use these lines as you wish. This is an efficient C implementation of the Ford-Fulkerson algorithm (Version of Edmonds-Karp) for maxflow. It should easily scale to hundreds of millions of edges if the data is not too adverserial.

# To compile:
"gcc maxflow.c -O3 -o maxflow".

# To execute:
"./maxflow edgelist.txt source target res.txt".
- "edgelist.txt" should contain the directed edges with capaities: one edge on each line (3 unsigned int separated by spaces "n1 n2 c").
- "source" and "target" are the source and target IDs.
- "res.txt" contains the results: "n1 n2 c f" on each line.

Will print some information in the terminal.


# Initial contributors

Maximilien Danisch  
May 2017  
http://bit.ly/maxdan94  
maximilien.danisch@gmail.com
