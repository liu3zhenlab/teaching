students <- c("Bilal Ahmad",
              "Ednaldo Borgato",
              "Heather Forster",
              "Tommy Galfano",
              "Candy Hernandez",
              "Molly Jones",
              "Afsana Noor",
              "Tyler Suelter",
              "Augusto Tessele",
              "Ben Wiens",
              "Xiaoting Xu")
nstudents <- length(students)
present_order <- sample(students, nstudents, replace=F)
data.frame(Order=1:nstudents, names=present_order )
