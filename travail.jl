# ---
# title: "Stratégie de vaccination spatiale : simulation d'une intervention épidémique sous contraintes"
# repository: ametrinee/BIO2045_devoir3
# auteurs:
#    - nom: Hong
#      prenom: Sumi
#      matricule: 20233311
#      github: ametrinee
#    - nom: Nivelle
#      prenom: Anna
#      matricule: 20349280
#      github: anna2808-prog
# ---

# # Introduction

# Les campagnes de vaccination à grande échelle ont permis de réduire la mortalité et d’atteindre, pour de nombreuses maladies infectieuses, 
# une immunité collective, c’est-à-dire la résistance d'une population à un agent contagieux lorsqu'une très grande proportion des individus 
# est immunisée contre celui-ci. En revanche, ces campagnes sont coûteuses en raison du prix d’une dose de vaccin et de l’administration d’un 
# test antigénique rapide (RAT) pour détecter les individus infectieux. Dans ce contexte, il apparaît pertinent de pouvoir prédire l’efficacité 
# d’une stratégie de vaccination, c’est-à-dire la réduction de la mortalité comparée à l’absence d’intervention, en fonction de son coût total.

# Le problème biologique central de cette étude est de contrôler la propagation d'une maladie infectieuse émergente
# au sein d'une population naïve, alors que le pathogène présente une dynamique de transmission
# silencieuse et une létalité absolue. L'enjeu est de concevoir une stratégie d'intervention épidémiologique capable
# de freiner la contagion et de minimiser la mortalité, tout en composant avec des informations incomplètes et des
# ressources matérielles et financières strictement limitées.

# Les modèles à base d’agent sont couramment utilisés dans le domaine de la prévention des épidémies car ils permettent de simuler chaque individu 
# séparément et donc de mieux intégrer les interactions sociales et les décisions individuelles @zhang2025agent.

# Dans ce travail, nous proposons un modèle à base d'agents qui intègre plusieurs contraintes biologiques et logistiques qui reflètent la réalité
# des crises sanitaires. La difficulté majeure réside dans la nature de la maladie : elle est strictement asymptomatique,
# se transmet par contact (taux de 0.4), et s'avère 100 % fatale à l'issue de 21 jours. Biologiquement, cela modélise un
# pathogène ayant une longue période de contagiosité silencieuse avant l'apparition de l'issue fatale. Comme l'ont
# démontré certains travaux @he2020temporal, la proportion de transmission pré-symptomatique ou asymptomatique est le facteur
# le plus déterminant qui rend une épidémie impossible à contrôler par de simples mesures d'isolement, nécessitant une 
# recherche active des cas. De plus, l'impossibilité d'intervenir avant le premier décès reflète les systèmes de surveillance
# basés sur la détection syndromique. Ce premier décès agit comme un signal d'alarme déclenchant la riposte sanitaire, à un
# stade où un réseau de transmission invisible s'est déjà établi.

# À cette dynamique s'ajoutent des contraintes matérielles sévères. Le dépistage s'effectue via des tests antigéniques
# rapides (RAT) présentant 95 % d'efficacité, ce taux d'erreur imposant une prise de décision sous incertitude. Par ailleurs,
# le vaccin nécessite deux jours pour conférer une protection, un délai biologiquement justifié par le temps nécessaire au
# système immunitaire adaptatif pour générer une réponse protectrice. L'ensemble de cette campagne doit être orchestré avec
# un budget strictement limité à 21 000 $, interdisant toute approche de vaccination de masse. Cette pénurie justifie le
# recours à une stratégie de dépistage localisé suivie d'une "vaccination en anneau" (Ring Vaccination), une approche
# historiquement validée pour briser les chaînes de transmission avec un nombre restreint de doses @henaorestrepo2017ebola.

# # Présentation du modèle

# Le modèle de base a été modifié pour intégrer les processus biologiques et logistiques de l'intervention sanitaire. D'un
# point de vue biologique, l'acquisition de l'immunité n'est pas instantanée. Pour refléter cette réalité, la structure
# représentant chaque individu a été enrichie de nouvelles variables pour suivre le statut vaccinal, le temps écoulé depuis
# la vaccination et l'historique de dépistage. Du point de vue algorithmique, la boucle de simulation a été mise à jour pour
# inclure une gestion stricte des ressources (décrémentation du budget à chaque test ou vaccin administré).

# Conformément aux contraintes, la campagne ne s'amorce qu'à l'enregistrement du premier décès. La stratégie adoptée repose
# sur le principe de la "vaccination en anneau" couplée à un dépistage ciblé. Concrètement, l'algorithme sélectionne aléatoirement
# des agents non testés et leur administre un test RAT. Lorsqu'un agent est identifié comme positif, une intervention spatiale
# immédiate est déclenchée : le modèle repère tous les individus présents dans la même cellule spatiale que le cas détecté.
# Si ces voisins ne sont pas encore vaccinés et que le budget le permet, ils reçoivent immédiatement une dose de vaccin. Cette
# méthode permet de construire de véritables barrières immunitaires autour des foyers infectieux pour étouffer la transmission spatiale.
 
# # Implémentation

import UUIDs
import Random
using StatsBase

# ## Définition des structures de données

## Structure représentant un individu et ses états (santé, vaccination, dépistage)
Base.@kwdef mutable struct Agent
    x::Int64 = 0
    y::Int64 = 0
    clock::Int64 = 21                # La durée de la maladie = 21 jours
    infectious::Bool = false
    id::UUIDs.UUID = UUIDs.uuid4()
    vaccinated::Bool = false         # Statut vaccinal
    days_after_vax::Int64 = 0        # Temps écoulé depuis la vaccination
    tested::Bool = false             # Historique de dépistage pour éviter les tests redondants
end

## Structure définissant les limites spatiales de la simulation
Base.@kwdef mutable struct Landscape
    xmin::Int64 = -50
    xmax::Int64 = 50
    ymin::Int64 = -50
    ymax::Int64 = 50
end

## Alias pour manipuler un groupe d'agents
const Population = Vector{Agent}

## Structure pour enregistrer les paramètres de chaque transmission
struct InfectionEvent
    time::Int64
    from::UUIDs.UUID
    to::UUIDs.UUID
    x::Int64
    y::Int64
end

# ## Fonctions du modèle

## Gère le mouvement aléatoire des agents sur la grille avec gestion des bordures
function move!(A::Agent, L::Landscape)
    A.x += rand(-1:1)
    A.y += rand(-1:1)
    ## Gestion des bordures
    A.x = clamp(A.x, L.xmin, L.xmax)
    A.y = clamp(A.y, L.ymin, L.ymax)
    return A
end

## Fonctions utilitaires pour filtrer la population selon leur état de santé
isinfectious(agent::Agent) = agent.infectious
ishealthy(agent::Agent) = !isinfectious(agent)
infectious(pop::Population) = filter(isinfectious, pop)
healthy(pop::Population) = filter(ishealthy, pop)

## Identifie les agents présents dans la même cellule (voisinage immédiat)
incell(target::Agent, pop::Population) = filter(ag -> (ag.x, ag.y) == (target.x, target.y), pop)

# ## Simulation principale

## Fonction cœur exécutant la dynamique épidémique et la stratégie d'intervention
function simulation(L::Landscape, n_initial::Int64, budget_total::Float64)
    pop = [Agent(x=rand(L.xmin:L.xmax), y=rand(L.ymin:L.ymax)) for _ in 1:n_initial]
    rand(pop).infectious = true
    
    tick = 0
    max_tick = 2000
    current_budget = budget_total
    events = InfectionEvent[]

    ## Listes pour le suivi temporel de la dynamique
    S_trace, I_trace, D_trace = Int64[], Int64[], Int64[]

    while (length(infectious(pop)) > 0) && (tick < max_tick)
        tick += 1
  
        ## Phase 1 : Déplacement des individus
        for agent in pop
            move!(agent, L)
        end
  
        ## Phase 2 : Intervention sanitaire (Dépistage RAT et vaccination)
        if length(pop) < n_initial && current_budget > 0
            for agent in Random.shuffle(pop)
                ## 4.0 $ = Coût du test RAT
                if current_budget >= 4.0 && !agent.tested
                    current_budget -= 4.0
                    agent.tested = true
     
                    ## Simulation du test RAT avec 95% d'efficacité
                    pos = agent.infectious ? rand() <= 0.95 : rand() <= 0.05

                    ## Si positif, vaccination des contacts dans la même cellule
                    if pos
                        for contact in incell(agent, pop)
                            ## 17.0 $ = Coût du vaccin
                            if current_budget >= 17.0 && !contact.vaccinated
                                current_budget -= 17.0
                                contact.vaccinated = true
                            end
                        end
                    end
                end
            end
        end
  
        ## Phase 3 : Transmission de la maladie
        for spreader in Random.shuffle(infectious(pop))
            for target in healthy(incell(spreader, pop))
                ## Vérification de l'immunité (active 2 jours après le vaccin)
                immune = target.vaccinated && target.days_after_vax >= 2
                ## 0.4 = Taux de transmission par contact
                if !immune && rand() <= 0.4
                    target.infectious = true
                    push!(events, InfectionEvent(tick, spreader.id, target.id, target.x, target.y))
                end
            end
        end
  
        ## Phase 4 : Évolution des états internes et mortalité
        for agent in pop
            if agent.infectious; agent.clock -= 1; end
            if agent.vaccinated; agent.days_after_vax += 1; end
        end
        ## Retrait des individus décédés
        pop = filter(a -> a.clock > 0, pop)
  
        ## Enregistrement des données de la génération
        push!(S_trace, length(healthy(pop)))
        push!(I_trace, length(infectious(pop)))
        push!(D_trace, n_initial - length(pop))
    end
    
    return S_trace, I_trace, D_trace, current_budget
end

# # Présentation des résultats

Random.seed!(123456)
using CairoMakie
using Statistics

# ## 1. Scénario de référence : Absence d'intervention

# Ce scénario sert de témoin négatif pour observer la dynamique naturelle du pathogène sans contraintes budgétaires ni mesures sanitaires.

L = Landscape()

## Exécution de la simulation avec un budget de 0 $
S_ref, I_ref, D_ref = simulation(L, 3750, 0.0)

## Création de la figure pour le scénario de référence
f_ref = Figure();
ax_ref = Axis(f_ref[1, 1], xlabel="Génération", ylabel="Population", title="Évolution de l'épidémie sans intervention")

## Tracé des courbes (S, I, D)
stairs!(ax_ref, S_ref, label="Sains", color=:midnightblue)
stairs!(ax_ref, I_ref, label="Infectieux", color=:red)
stairs!(ax_ref, D_ref, label="Décédés", color=:dimgrey)

## Légende
axislegend(ax_ref)

save("sans-intervention.png", f_ref)
f_ref

# **Figure 1. Dynamique épidémique de référence sans intervention sanitaire.**

# ## 2. Simulation de la stratégie : Intervention avec vaccination

# Ce scénario évalue l'efficacité de la vaccination en anneau sous une contrainte budgétaire stricte de 21 000 $.

L = Landscape()

## Exécution de la simulation avec le budget alloué (21 000$)
S, I, D, final_budget = simulation(L, 3750, 21000.0)

## Génération du graphique de la dynamique épidémique sous intervention
f = Figure();
ax = Axis(f[1, 1], xlabel="Génération", ylabel="Population", title="Dynamique de l'épidémie avec vaccination")

## Tracé des courbes (S, I, D)
stairs!(ax, S, label="Sains", color=:midnightblue)
stairs!(ax, I, label="Infectieux", color=:red)
stairs!(ax, D, label="Décédés", color=:dimgrey)

## Légende
axislegend(ax)

save("avec-intervention.png", f)
f

# **Figure 2. Impact de la vaccination en anneau sur la stabilisation de la population saine.**

# La comparaison entre les deux figures montre une nette réduction de la mortalité grâce à l'intervention.
# En l'absence de mesures, la population saine s'effondre rapidement. Avec notre stratégie, la courbe des
# survivants se stabilise à un niveau significativement plus élevé, prouvant que la vaccination en anneau
# parvient à isoler les foyers infectieux malgré un budget limité.

# # Réplications et analyse de la variabilité

# ## 1. Initialisation des conteneurs de données
n_reps = 100
results_S = Int64[];         # Population saine finale
results_D = Int64[];         # Nombre de décès final
results_budget = Float64[];  # Budget restant

# ## 2. Exécution de la boucle de simulation

Random.seed!(123456)

for i in 1:n_reps
    ## Exécution répétée pour capturer la variabilité stochastique du modèle
    S_trace, I_trace, D_trace, final_budget_tmp = simulation(L, 3750, 21000.0)

    ## Stockage des valeurs finales de la dernière génération dans les conteneurs
    push!(results_S, S_trace[end])
    push!(results_D, D_trace[end])
    push!(results_budget, final_budget_tmp)
end

# ## 3. Calcul des statistiques descriptives

mean_S = mean(results_S);             # Estimation de l'espérance du nombre de survivants
std_S = std(results_S);               # Mesure de la variabilité des survivants
mean_D = mean(results_D);             # Nombre moyen de décès sur 100 réplications
std_D = std(results_D);               # Variabilité de la mortalité entre les simulations
mean_budget = mean(results_budget);   # Budget moyen résiduel après l'intervention
std_budget = std(results_budget);     # Écart-type du coût financier de la campagne

# ## 4. Affichage des résultats pour l'analyse de robustesse

println("Moyenne des survivants : $mean_S (± $std_S)")
println("Moyenne des décès : $mean_D (± $std_D)")
println("Budget restant moyen : $mean_budget \$ (± $std_budget \$)")

# ## 5. Visualisation simple du budget

f_dots = Figure();
ax_dots = Axis(f_dots[1, 1], title = "Budget restant (100 réps)", ylabel = "Budget (\$)")
scatter!(ax_dots, results_budget, markersize = 12, color = (:midnightblue, 0.5))
save("budget-dots.png", f_dots)
f_dots

# **Figure 3. Analyse de la variabilité stochastique : distribution bimodale du budget résiduel.**

# La Figure ci-dessus illustre la répartition des résultats financiers.
# On observe une distribution bimodale : dans la grande majorité des cas (92 %),
# le budget est presque entièrement consommé pour freiner l'épidémie,
# tandis que dans une minorité de simulations (8 %), il reste proche
# de 21 000 $. Cette répartition stochastique se traduit par une moyenne
# de 3 101 $ avec un écart-type de 5 317 $.

# # Discussion

# Cette étude vise à évaluer l’efficacité d’une stratégie de vaccination en anneau sous contrainte budgétaire.
# Pour ce faire, un modèle à base d’agents a été utilisé pour simuler la dynamique épidémique.

# Les résultats montrent une éradication complète de l’épidémie avant l’extinction de la population. La vaccination
# en anneau a considérablement réduit la mortalité et permis d’atteindre une immunité locale permettant de rompre
# les chaînes de transmission. Cependant, le nombre de survivants présente une variabilité importante, avec un
# écart-type de 805, probablement en raison de la stochasticité du modèle, en particulier la dynamique spatiale
# aléatoire. Le modèle pourrait être plus réaliste en créant une population hétérogène et des coordonnées spécifiques
# plus susceptibles d’être fréquentées pour imiter des zones de haute densité telles que les lieux de travail, les
# écoles ou les magasins, éventuellement avec des règles de distanciation physique. De plus, il pourrait être pertinent
# de réduire la probabilité d’interaction avec des individus infectés et testés pour prendre en compte les changements
# de comportement des individus lorsqu’ils sont entourés de malades, ou bien les confinements pouvant être imposés face
# à certaines épidémies comme celle de la Covid-19 en 2020 @he2020temporal. Par ailleurs, l'impossibilité d'intervenir
# avant le premier décès crée une fenêtre de transmission invisible. Cette latence administrative est un facteur critique
# qui explique pourquoi une partie de la population succombe avant l'établissement des barrières immunitaires.

# Le modèle est également simplifié d’un point de vue biologique. Une durée de 21 jours de la maladie a été considérée,
# bien que celle-ci puisse en réalité varier en fonction de l'âge, la génétique ou encore la charge virale @langford2007predictors.
# L’hypothèse selon laquelle l’infection est toujours létale pourrait également être reconsidérée selon le micro-organisme
# étudié. Il est aussi à noter que la protection totale du vaccin après 2 jours est un cas idéal, la réponse vaccinale ne
# pouvant pas être sûre à 100 % et diminuant avec le temps en raison de la perte d’anticorps @dalmat2021efficacite. Un autre
# aspect négligé est la période d’incubation, durant laquelle les individus peuvent être infectés mais pas contagieux.

# Concernant le coût de la campagne, la stratégie de vaccination en anneau a consommé en moyenne 17 899 $, soit une économie
# massive par rapport à une vaccination de masse (3 750 × 21 $ = 78 750 $). Cet écart-type élevé de 5 317 $ s'explique par la
# distribution bimodale. Les résultats se divisent en deux trajectoires opposées : dans 8 % des cas, l'épidémie s'éteint par
# simple stochasticité dès le début, préservant ainsi le budget, tandis que dans les 92 % restants, le virus se propage davantage,
# forçant une intervention massive qui consomme la quasi-totalité des fonds. Cette absence de scénario moyen explique la forte
# variabilité des données, tant pour les coûts que pour la mortalité.

# Malgré ces limites, ce modèle met en évidence la pertinence de la vaccination en anneau lorsque les ressources sont limitées
# face à des épidémies asymptomatiques @henaorestrepo2017ebola.

