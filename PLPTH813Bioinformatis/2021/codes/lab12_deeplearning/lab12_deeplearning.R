# A simple deep learning example
# Sanzhen Liu
# PLPTH813
# 4/29/2021
# materials were from https://tensorflow.rstudio.com/tutorials/beginners/basic-ml/tutorial_basic_classification/
library(tensorflow)
library(keras)
fashion_mnist <- dataset_fashion_mnist()
summary(fashion_mnist$train)
c(train_images, train_labels) %<-% fashion_mnist$train
c(test_images, test_labels) %<-% fashion_mnist$test

class_names = c('T-shirt/top',
                'Trouser',
                'Pullover',
                'Dress',
                'Coat',
                'Sandal',
                'Shirt',
                'Sneaker',
                'Bag',
                'Ankle boot')

dim(train_images)
dim(train_labels)

library(tidyr)
library(ggplot2)

image_1 <- as.data.frame(train_images[1, , ])
colnames(image_1) <- seq_len(ncol(image_1))
image_1$y <- seq_len(nrow(image_1))
image_1 <- gather(image_1, "x", "value", -y)
image_1$x <- as.integer(image_1$x)

ggplot(image_1, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", na.value = NA) +
  scale_y_reverse() +
  theme_minimal() +
  theme(panel.grid = element_blank())   +
  theme(aspect.ratio = 1) +
  xlab("") +
  ylab("")

#We scale these values to a range of 0 to 1 before feeding to the neural network model. For this, we simply divide by 255.
#It’s important that the training set and the testing set are preprocessed in the same way:
train_images <- train_images / 255
test_images <- test_images / 255

## display top 9 images
par(mfcol=c(3,3))
par(mar=c(0, 0, 1.5, 0), xaxs='i', yaxs='i')
for (i in 1:9) {
  img <- train_images[i, , ]
  img <- t(apply(img, 2, rev))
  image(1:28, 1:28, img, col = gray((0:255)/255), xaxt = 'n', yaxt = 'n',
        main = paste(class_names[train_labels[i] + 1]))
}

### setup layers
library(keras)
model <- keras_model_sequential() %>%
  layer_flatten(input_shape = c(28, 28)) %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dense(units = 10, activation = 'softmax')

### compile
model %>% compile(
  optimizer = 'adam',
  loss = 'sparse_categorical_crossentropy',
  metrics = c('accuracy')
)

### training
model %>% fit(train_images, train_labels, epochs=20,
              validation_split=0.3, verbose=2)

### select epoch=5
model %>% fit(train_images, train_labels, epochs=5, verbose=2)

### prediction
predictions <- model %>% predict(test_images)
predictions[1, ]
which.max(predictions[1, ])
class_names[which.max(predictions[1, ])]

### manual check
par(mfcol=c(4,4))
par(mar=c(0, 0, 1.5, 0), xaxs='i', yaxs='i')
for (i in 1:16) {
  img <- test_images[i, , ]
  img <- t(apply(img, 2, rev))
  # subtract 1 as labels go from 0 to 9
  predicted_label <- which.max(predictions[i, ]) - 1
  true_label <- test_labels[i]
  if (predicted_label == true_label) {
    color <- '#008800'
  } else {
    color <- '#bb0000'
  }
  image(1:28, 1:28, img, col = gray((0:255)/255), xaxt = 'n', yaxt = 'n',
        main = paste0(class_names[predicted_label + 1], " (",
                      class_names[true_label + 1], ")"),
        col.main = color)
}

### predict one image
# Grab an image from the test dataset
img <- test_images[1, , , drop = FALSE]
dim(img)
prediction1 <- model %>% predict(img)
prediction1
