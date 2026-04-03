# ---
# title: "Stratégie de vaccination spatiale : simulation d'une intervention épidémique sous contraintes"
# repository: ametrinee/BIO2045_devoir3
# auteurs:
#    - nom: Auteur
#      prenom: Premier
#      matricule: XXXXXXXX
#      github: premierAuteur
#    - nom: Auteur
#      prenom: Deuxième
#      matricule: XXXXXXXX
#      github: DeuxiAut
# ---

# # Introduction

# La gestion des épidémies représente l'un des défis majeurs de la santé publique et de l'épidémiologie spatiale.
# Le problème biologique central de cette étude est de contrôler la propagation d'une maladie infectieuse émergente
# au sein d'une population naïve (3750 individus), alors que le pathogène présente une dynamique de transmission
# silencieuse et une létalité absolue. L'enjeu est de concevoir une stratégie d'intervention épidémiologique capable
# de freiner la contagion et de minimiser la mortalité, tout en composant avec des informations incomplètes et des
# ressources matérielles et financières strictement limitées.

# Pour simuler ce défi, notre modèle intègre plusieurs contraintes biologiques et logistiques qui reflètent la réalité
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

# Structure représentant un individu et ses états (santé, vaccination, dépistage)
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

# Structure définissant les limites spatiales de la simulation
Base.@kwdef mutable struct Landscape
    xmin::Int64 = -50
    xmax::Int64 = 50
    ymin::Int64 = -50
    ymax::Int64 = 50
end

# Alias pour manipuler un groupe d'agents
const Population = Vector{Agent}

# Structure pour enregistrer les paramètres de chaque transmission
struct InfectionEvent
    time::Int64
    from::UUIDs.UUID
    to::UUIDs.UUID
    x::Int64
    y::Int64
end

# ## Fonctions du modèle

# Gère le mouvement aléatoire des agents sur la grille avec gestion des bordures
function move!(A::Agent, L::Landscape)
    A.x += rand(-1:1)
    A.y += rand(-1:1)
    # Gestion des bordures
    A.x = clamp(A.x, L.xmin, L.xmax)
    A.y = clamp(A.y, L.ymin, L.ymax)
    return A
end

# Fonctions utilitaires pour filtrer la population selon leur état de santé
isinfectious(agent::Agent) = agent.infectious
ishealthy(agent::Agent) = !isinfectious(agent)
infectious(pop::Population) = filter(isinfectious, pop)
healthy(pop::Population) = filter(ishealthy, pop)

# Identifie les agents présents dans la même cellule (voisinage immédiat)
incell(target::Agent, pop::Population) = filter(ag -> (ag.x, ag.y) == (target.x, target.y), pop)

# ## Simulation principale

# Fonction cœur exécutant la dynamique épidémique et la stratégie d'intervention
function run_epidemic_simulation(L::Landscape, n_initial::Int64, budget_total::Float64)
    pop = [Agent(x=rand(L.xmin:L.xmax), y=rand(L.ymin:L.ymax)) for _ in 1:n_initial]
    rand(pop).infectious = true
    
    tick = 0
    max_tick = 2000
    current_budget = budget_total
    events = InfectionEvent[]

    # Listes pour le suivi temporel de la dynamique
    S_trace, I_trace, D_trace = Int64[], Int64[], Int64[]

    while (length(infectious(pop)) > 0) && (tick < max_tick)
        tick += 1
  
        # Phase 1 : Déplacement des individus
        for agent in pop
            move!(agent, L)
        end
  
        # Phase 2 : Intervention sanitaire (Dépistage RAT et vaccination)
        if length(pop) < n_initial && current_budget > 0
            for agent in Random.shuffle(pop)
                if current_budget >= 4.0 && !agent.tested
                    current_budget -= 4.0
                    agent.tested = true
     
                    # Simulation du test RAT avec 95% d'efficacité
                    pos = agent.infectious ? rand() <= 0.95 : rand() <= 0.05

                    # Si positif, vaccination des contacts dans la même cellule
                    if pos
                        for contact in incell(agent, pop)
                            if current_budget >= 17.0 && !contact.vaccinated
                                current_budget -= 17.0
                                contact.vaccinated = true
                            end
                        end
                    end
                end
                current_budget < 4.0 && break
            end
        end
  
        # Phase 3 : Transmission de la maladie
        for spreader in Random.shuffle(infectious(pop))
            for target in healthy(incell(spreader, pop))
                # Vérification de l'immunité (active 2 jours après le vaccin)
                immune = target.vaccinated && target.days_after_vax >= 2
                if !immune && rand() <= 0.4
                    target.infectious = true
                    push!(events, InfectionEvent(tick, spreader.id, target.id, target.x, target.y))
                end
            end
        end
  
        # Phase 4 : Évolution des états internes et mortalité
        for agent in pop
            if agent.infectious; agent.clock -= 1; end
            if agent.vaccinated; agent.days_after_vax += 1; end
        end
        # Retrait des individus décédés
        pop = filter(a -> a.clock > 0, pop)
  
        # Enregistrement des données de la génération
        push!(S_trace, length(healthy(pop)))
        push!(I_trace, length(infectious(pop)))
        push!(D_trace, n_initial - length(pop))
    end
    
    return S_trace, I_trace, D_trace, current_budget
end

# # Présentation des résultats

import Random
Random.seed!(123456)
using CairoMakie
using Statistics

# ## Simulation unique et visualisation
L = Landscape()
S, I, D, final_budget = run_epidemic_simulation(L, 3750, 21000.0)

# Génération du graphique de la dynamique épidémique
f = Figure();
ax = Axis(f[1, 1], xlabel="Génération", ylabel="Population", title="Dynamique de l'épidémie avec vaccination")
stairs!(ax, S, label="Sains", color=:black)
stairs!(ax, I, label="Infectieux", color=:red)
stairs!(ax, D, label="Décédés", color=:grey)
axislegend(ax)
save("travail-simulation.png", f)
f

# ## Réplications et analyse de la variabilité

# 1. Configuration du nombre de réplications
n_reps = 10
results_S = Int64[]    # Population saine finale
results_D = Int64[]    # Nombre de décès final

for i in 1:n_reps
    # Exécution répétée pour capturer la variabilité stochastique
    S_trace, I_trace, D_trace, final_budget = run_epidemic_simulation(L, 3750, 21000.0)

    # Stockage uniquement des valeurs finales de la dernière génération
    push!(results_S, S_trace[end])
    push!(results_D, D_trace[end])
end

# 2. Calcul des statistiques globales
mean_S = mean(results_S)
std_S = std(results_S)
mean_D = mean(results_D)

println("Moyenne des survivants : $mean_S (Écart-type : $std_S)")

# # Discussion
