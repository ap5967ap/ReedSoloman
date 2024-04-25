
# Reed Soloman

This is a repository containing the implementation of integer version of Reed Soloman Error correcting code.
The code works for input message upto 10¹⁰⁰⁰. I have used gmpy2 for implementation of big integers.





## Run Locally



```bash
 git clone https://github.com/ap5967ap/ReedSoloman.git
 cd ReedSoloman
```

Installing dependencies
```bash
 pip install gmpy2
```

Running 
```bash
 python3 reed_soloman.py
```    
## Running Tests

To run for multiple rounds to give accuracy details

```bash
  python3 tst.py
```

The code returns the reconstructed message if it reconstructs the message correctly else it returns -1;

## Authors

- [@ap5967ap](https://www.github.com/ap5967ap)

