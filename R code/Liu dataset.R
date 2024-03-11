

########handle Liu dataset
counts = read.csv("./data/s_data/Liu/Liu_tpm.txt", sep = "\t")
patients_name = counts[,1]
counts= counts[,-1]
counts = as.data.frame(t(counts))
colnames(counts) = patients_name
tpm = cbind(data.frame(gene_id = rownames(counts)), counts)
write.table(tpm, file = "./data/h_data/Liu/Liu_tpm.txt", quote = F, sep = "\t", row.names = F)

tpm = log2(tpm + 1)
tpm = cbind(data.frame(gene_id = rownames(tpm)), tpm)
write.table(tpm, file = "./data/h_data/Liu/Liu_log2tpm.txt", quote = F, sep = "\t", row.names = F)


########handle the labels
response_labels = read.csv("./data/s_data/Liu/Liu_patient.txt", sep = "\t", header = F)
patients_name = response_labels$V1
flags = vector()
response = vector()
for (i in 1:length(response_labels[,1]))
{
  if (response_labels[i, 2] == "responder")
  {
    flags = append(flags, 1)
    response = append(response, "R")
  }
  else
  {
    flags = append(flags, 0)
    response = append(response, "NR")
  }
}
response_labels = data.frame(Patient = patients_name, Response = response, flag = flags)
write.table(response_labels, file = "./data/h_data/Liu/Liu_patient.txt", quote = F, sep = "\t", row.names = F)
