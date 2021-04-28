### install tensorflow in R
install.packages("tensorflow")
library("tensorflow")

install_tensorflow(
  method = "conda",
  conda_python_version = "3.6",
)

### confirm that the installation succeeded with:
library("tensorflow")
tf$constant("Hellow Tensorflow")

### install keras
install.packages("keras")
