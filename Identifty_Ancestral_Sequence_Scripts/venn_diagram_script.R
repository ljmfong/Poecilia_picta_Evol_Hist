#Shared Kmers Venn Diagram

library(VennDiagram)

#From the numbers taken from HAWK

### With Parae:
#area 1 = P. ret
#area 2 = P. win
#area 3 = P. picta
#area 4 = P. parae
#draw.quad.venn(area1, area2, area3, area4, 
#n12, n13, n14, n23, n24,n34, 
#n123, n124, n134, n234, n1234, category = rep("",4), lwd = rep(2, 4), 
#lty = rep("solid", 4), col = rep("black", 4)

#MALE UNIQUE:

draw.quad.venn(1220830,1122342,1137779,960495,
               20947,1071,771,243,129,119787,
               4,21,430,0,0,
               category = rep("", 4),
               lwd = rep(1,4), lty = rep("solid", 4), 
               col = rep("midnightblue", 4), fill = "dodgerblue",
               alpha = rep(0.2, 4), label.col = rep("black", 15), cex = rep(0.8, 15),
               fontface = rep("plain", 15), fontfamily = rep("serif", 15))

#draw.quad.venn(area1, area2, area3, area4, 
#n12, n13, n14, n23, n24,n34, 
#n123, n124, n134, n234, n1234)

#FEMALE UNIQUE:

draw.quad.venn(958631,878688,545336,503322,
               484,21,85,104,10,325,
               0,0,0,0,0,
               category = rep("", 4),
               lwd = rep(1,4), lty = rep("solid", 4), 
               col = rep("maroon", 4), fill = "red",
               alpha = rep(0.2, 4), label.col = rep("black", 15), cex = rep(0.8, 15),
               fontface = rep("plain", 15), fontfamily = rep("serif", 15))





