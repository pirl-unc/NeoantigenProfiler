library(lpSolve)

set.seed(1)

sim_antigens <- function(n, m, k) {
	tab1 <- lapply(1:n, function(j) {
		ifelse(1:m %in% sample(1:m, k), 1, 0)
	})
	tab2 <- lapply(1:m, function(j) {
		res <- rep(0, m)
		res[j] <- 1
		res
	})
	c(tab1, tab2, tab2)	
}

run <- function(dat, weights) {
	f.obj <- weights
	f.con <- dat
	len1 <- nrow(dat) - 2 * length(weights)
	len2 <- length(weights)
	f.dir <- c(rep('>=', len1 + len2), rep('<=', len2))
	f.rhs <- c(rep(1, len1), rep(0, len2), rep(1, len2))

	solution <- lp ("min", f.obj, f.con, f.dir, f.rhs)
	res <- solution$solution
	npick <- runif(len2)
	pick <- ifelse(res >= npick, 1, 0)
	pick
}

run1 <- function(dat, weights) {
	f.obj <- weights
	len1 <- nrow(dat) - 2 * length(weights)
	len2 <- length(weights)
	f.con <- dat[1:len1, 1:len2]
	f.dir <- c(rep('>=', len1))
	f.rhs <- c(rep(1, len1))

	solution <- lp ("min", f.obj, f.con, f.dir, f.rhs)
	res <- solution$solution
	npick <- runif(len2)
	pick <- ifelse(res >= npick, 1, 0)
	pick
}

iterate <- function(dat, weights, times) {
	len1 <- nrow(dat) - 2 * length(weights)
	len2 <- length(weights)
	runs <- lapply(1:times, function(i) {
		run(dat, weights)
	})
	runs <- do.call(rbind, runs)
	runs <- colSums(runs)
	res <- ifelse(runs > 0, 1, 0)
	dat100 <- dat[1:len1, 1:len2]
	message("Selected ", sum(res), " antigens.")
	message("Total Covered: ", sum(rowSums(dat100[,res==1])>0))
	message("Total Weight: ", sum(weights * res))
	list(which(res>0), which(rowSums(dat100[,res==1]) > 0))
}

run_sim <- function(ncells, nantigens, k, times) {
	dat <- sim_antigens(ncells, nantigens, k)
	dat <- do.call(rbind, dat)
	weights <- runif(nantigens)
	iterate(dat, weights, times)
}

#run_sim(1000, 300, 10, 3)

ncells <- 1000
nantigens <- 100
k <- 30

dat <- sim_antigens(ncells, nantigens, k)
dat <- do.call(rbind, dat)
weights <- runif(nantigens)

iter <- 100

res <- lapply (1:iter, function(i) {

	coverage <- c()
	antigens <- c()
	times <- 0

	while(length(coverage) != ncells) {
		message("Start coverage: ", length(coverage))
		li <- iterate(dat, weights, 1)
		antigens <- c(antigens, as.character(li[[1]]))
		antigens <- antigens[!duplicated(antigens)]
		coverage <- c(coverage, as.character(li[[2]]))
		coverage <- coverage[!duplicated(coverage)]
		times <- times + 1
		message("Coverage: ", length(coverage))
		message("Antigens: ", length(antigens))
		message("Times: ", times)
	}
	list(coverage, antigens, times)
})

coverages <- vapply(1:iter, function(i) {length(res[[i]][[1]])}, numeric(1))
antigens <- vapply(1:iter, function(i) {length(res[[i]][[2]])}, numeric(1))
times <- vapply(1:iter, function(i) {res[[i]][[3]]}, numeric(1))

ha <- hist(antigens, breaks=(max(antigens) - min(antigens)))
ht <- hist(times, breaks=(max(times) - min(times)))


antis <- lapply(1:iter, function(i) {res[[i]][[2]]})

tab <- table(unlist(antis))
tab <- tab[order(tab, decreasing=T)]

anti_df <- data.frame(tab)
anti_df[,1] <- as.character(anti_df[,1])
colnames(anti_df) <- c("Antigen", "Occurences")

g <- ggplot(data=anti_df, aes(x=Antigen, y=Occurences)) +
	geom_bar(stat="identity", aes(x=factor(Antigen, level=anti_df[,1]))) +
	ggtitle("Frequnecy of selection of Antigens")

#iterate(dat, weights, times)
