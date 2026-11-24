# Maquette 3D du centre de Caen

Maquette 3D interactive du centre-ville de Caen (2,4 × 1,8 km), dans le navigateur, avec Three.js.

## Contenu

- `caen/index.html` : la page (Three.js 0.147, aucune étape de build).
- `caen/caen.json` : les données de la maquette, générées depuis OpenStreetMap.
- `data/osm.json`, `data/extra.json` : extraits Overpass bruts (bâti, voirie, eau, parcs, arbres, tram, mobilier).
- `data/tram_platforms.json` : quais du tram.
- `scripts/convert.py` : conversion des extraits OSM en `caen.json`.

## Fonctions

- Heure du jour avec soleil, crépuscule et nuit (fenêtres, lampadaires et monuments éclairés).
- Toits d'ardoise à brisis avec lucarnes et cheminées, façades calées étage par étage (Reconstruction, vieux Caen, pans de bois, moderne).
- Château, églises et monuments civils détaillés ; mode « jeu vidéo » low-poly avec nuit.
- Eau miroir, pluie, nuages, piétons, cyclistes, circulation, terrasses.
- Visite guidée, carte postale, curseur d'époque (été 1944 → aujourd'hui, reconstitution indicative).
- Quartier du Vaugueux reconstitué maison par maison, presqu'île détaillée bâtiment par bâtiment.
- Tramway : voies, ligne aérienne et vraies stations ; eau creusée sous les quais, ponts en 3D.

## Lancer en local

```
python -m http.server 8765 -d caen
```

puis ouvrir http://localhost:8765.

## Régénérer les données

```
python scripts/convert.py data/osm.json data/extra.json caen/caen.json
```

## Données

- Bâti, voirie, eau, végétation, tram : © contributeurs OpenStreetMap, sous licence ODbL.
