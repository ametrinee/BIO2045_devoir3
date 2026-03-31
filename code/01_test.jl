# ### Test

# Le contenu de ce fichier est dans `code/01_test.jl` -- il est inclus dans le
# fichier principal, et dans le document PDF.

# Vous pouvez utiliser autant de fichiers dans le dossier `code` que nécessaire,
# et vous pouvez les nommer comme vous voulez. Il faut simplement vous assurer
# de les inclure dans le fichier `travail.jl`.

# Les figures qui sont créés dans un fichier inclus seront automatiquement
# ajoutées au PDF.

f = Figure()
ax = Axis(f[1,1]; aspect=DataAspect())
scatter!(ax, sort(randn(10)), sort(randn(10)), color=:black)
current_figure()

# Vous pouvez, par exemple, créer une fichier pour chaque figure des résultats.