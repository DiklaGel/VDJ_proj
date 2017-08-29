data = read.csv("final_table.csv")
png('hist_clones.png', width = 1600, height = 1200)
ggplot(d,aes(x=clone_id,group=Patient,fill=Patient))+ geom_histogram(position="identity",bins = length(levels(factor(d$clone_id))))+theme_bw()
dev.off()
png('hist_clones_per_patient.png', width = 1704, height = 2272)
ggplot(d,aes(x=clone_id, color=Patient))+geom_histogram(bins = length(levels(factor(d$clone_id))))+facet_wrap(~ Patient,nrow = 3)+theme_bw()
dev.off()
png('most_abundant_clones.png', width = 1600, height = 1200)
ggplot(d,aes(x=clone_id,group=Patient,fill=Patient))+ geom_histogram(position="identity",bins = length(levels(factor(d$clone_id))))+theme_bw()
dev.off()

vector = apply(matrix(1:max(data$clone_id)),1,function(x) sum(data[which(data$clone_id == x),"reads"]))

for(i in 1:length(data$clone_id)){
  data$size_realtive_to_clone[i] = data$reads[i]/vector[data$clone_id[i]]
}

data$size_realtive_to_clone = apply(matrix(1:length(data$clone_id)), 1, function(i) data$reads[i]/vector[data$clone_id[i]])

sum = 0
for(i in 2:length(out)){
  if(length(unique(out[[i]]$Patient))>1){
    print(paste("clone id: ", i))
    print(out[[i]])
    sum = sum + 1
  }
}
print(sum)



