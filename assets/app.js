const STRING_API = "https://string-db.org/api/tsv/interaction_partners";
const STRING_LINK = "https://string-db.org/network/";
const STRING_LIMIT = 1000;
const CHRONOS_STEP = 0.05;
const DEFAULT_CHRONOS_BOUNDS = { max: 1, min: -2 };
const SVG_NS = "http://www.w3.org/2000/svg";
const evidenceLabels = {
  nscore: "neighborhood",
  fscore: "fusion",
  pscore: "cooccur.",
  ascore: "coexpression",
  escore: "experiments",
  dscore: "databases",
  tscore: "textmining",
};

const state = {
  cache: null,
  center: "TP53",
  chronos: new Map(),
  chronosBounds: { ...DEFAULT_CHRONOS_BOUNDS },
  chronosMax: DEFAULT_CHRONOS_BOUNDS.max,
  dependencyScoreThreshold: -1.0,
  edges: [],
  essential: new Set(),
  ipdbEvidence: new Map(),
  ipdbPairs: new Map(),
  ipdbProteins: new Map(),
  knownLigands: new Map(),
  knownLigandOnly: false,
  physicalOnly: false,
  scoreMin: 0.99,
  scoreMax: 1,
  onlyEssential: true,
  activeCard: null,
  selected: null,
  source: "STRING API",
};

const els = {
  form: document.querySelector("#search-form"),
  input: document.querySelector("#protein-input"),
  suggestions: document.querySelector("#protein-suggestions"),
  graphWrap: document.querySelector(".graph-wrap"),
  graph: document.querySelector("#ppi-graph"),
  title: document.querySelector("#graph-title"),
  empty: document.querySelector("#empty-state"),
  knownLigandOnly: document.querySelector("#known-ligand-only"),
  loading: document.querySelector("#loading-pill"),
  apply: document.querySelector("#apply-filters"),
  min: document.querySelector("#min-score"),
  max: document.querySelector("#max-score"),
  minValue: document.querySelector("#min-score-value"),
  maxValue: document.querySelector("#max-score-value"),
  chronosMax: document.querySelector("#chronos-max"),
  chronosMaxValue: document.querySelector("#chronos-max-value"),
  essentialOnly: document.querySelector("#essential-only"),
  physicalOnly: document.querySelector("#physical-only"),
  targetCard: document.querySelector("#target-card"),
  targetClose: document.querySelector("#target-card-close"),
  targetName: document.querySelector("#target-name"),
  targetSections: document.querySelector("#target-card-sections"),
  visibleCount: document.querySelector("#visible-count"),
};

init();

async function init() {
  await Promise.all([
    loadEssentialProteins(),
    loadCache(),
    loadIpdbData(),
    loadKnownLigands(),
  ]);
  configureChronosControls();
  populateSuggestions();

  els.form.addEventListener("submit", (event) => {
    event.preventDefault();
    searchCurrentProtein();
  });
  els.apply.addEventListener("click", searchCurrentProtein);
  els.essentialOnly.addEventListener("change", () => {
    state.onlyEssential = els.essentialOnly.checked;
    state.selected = null;
    render();
  });
  els.knownLigandOnly.addEventListener("change", () => {
    state.knownLigandOnly = els.knownLigandOnly.checked;
    state.selected = null;
    render();
  });
  els.physicalOnly.addEventListener("change", () => {
    state.physicalOnly = els.physicalOnly.checked;
    state.selected = null;
    searchCurrentProtein();
  });
  els.min.addEventListener("input", syncScoreControls);
  els.max.addEventListener("input", syncScoreControls);
  els.chronosMax.addEventListener("input", syncChronosControls);
  els.targetClose.addEventListener("click", closeProteinCard);
  els.graph.addEventListener("click", selectNode);
  els.graph.addEventListener("keydown", selectNode);
  window.addEventListener("resize", render);

  await searchCurrentProtein();
}

async function loadEssentialProteins() {
  try {
    const response = await fetch("data/essential_proteins.json", {
      cache: "no-store",
    });
    const payload = await response.json();
    state.essential = new Set((payload.proteins || []).map(normalizeProtein));
    state.chronos = new Map(
      Object.entries(payload.summary || {}).map(([symbol, summary]) => [
        normalizeProtein(symbol),
        summary,
      ]),
    );
    state.chronosBounds = chronosBoundsFromPayload(payload);
    state.dependencyScoreThreshold = payload.meta?.score_threshold ?? -1.0;
  } catch (error) {
    console.warn("DepMap dependency dataset was not loaded", error);
    state.chronos = new Map();
    state.essential = new Set();
    state.chronosBounds = { ...DEFAULT_CHRONOS_BOUNDS };
  }
}

async function loadCache() {
  try {
    const response = await fetch("data/string_interactions.json", {
      cache: "no-store",
    });
    state.cache = await response.json();
  } catch (error) {
    console.warn("Local STRING cache was not loaded", error);
    state.cache = { interactions: [] };
  }
}

async function loadIpdbData() {
  try {
    const response = await fetch("data/ipdb_processed.json", {
      cache: "no-store",
    });
    const payload = await response.json();
    state.ipdbEvidence = normalizeEntryMap(payload.evidence || {});
    state.ipdbPairs = new Map(
      Object.entries(payload.pairs || {}).map(([key, pair]) => {
        const [left, right] = key.split("::");
        return [normalizePairKey(left, right), pair];
      }),
    );
    state.ipdbProteins = normalizeEntryMap(payload.proteins || {});
  } catch (error) {
    console.warn("Processed IPDB data was not loaded", error);
    state.ipdbEvidence = new Map();
    state.ipdbPairs = new Map();
    state.ipdbProteins = new Map();
  }
}

async function loadKnownLigands() {
  try {
    const response = await fetch("data/known_ligands.json", {
      cache: "no-store",
    });
    const payload = await response.json();
    state.knownLigands = new Map(
      Object.entries(payload.ligands || {}).map(([symbol, binding]) => [
        normalizeProtein(symbol),
        binding,
      ]),
    );
  } catch (error) {
    console.warn("BindingDB known ligand index was not loaded", error);
    state.knownLigands = new Map();
  }
}

function populateSuggestions() {
  const names = new Set([...state.chronos.keys(), ...state.essential]);

  state.cache.interactions.forEach((edge) => {
    names.add(edge.preferredName_A);
    names.add(edge.preferredName_B);
  });

  els.suggestions.replaceChildren(
    ...[...names]
      .filter(Boolean)
      .sort()
      .slice(0, 120)
      .map((name) => {
        const option = document.createElement("option");
        option.value = name;
        return option;
      }),
  );
}

function configureChronosControls() {
  const { max, min } = state.chronosBounds;

  els.chronosMax.disabled = state.chronos.size === 0;
  els.chronosMax.min = min.toFixed(2);
  els.chronosMax.max = max.toFixed(2);
  els.chronosMax.step = CHRONOS_STEP.toFixed(2);

  const defaultMax = clamp(state.dependencyScoreThreshold, min, max);
  els.chronosMax.value = defaultMax.toFixed(2);
  state.chronosMax = defaultMax;
  renderChronosControls();
}

function chronosBoundsFromPayload(payload) {
  const meanScores = Object.values(payload.summary || {})
    .map((summary) => Number(summary.mean_score))
    .filter(Number.isFinite);
  const metaMin = Number(payload.meta?.mean_score_range?.min);
  const metaMax = Number(payload.meta?.mean_score_range?.max);

  if (!meanScores.length && !Number.isFinite(metaMin) && !Number.isFinite(metaMax)) {
    return { ...DEFAULT_CHRONOS_BOUNDS };
  }

  const dataMin = Number.isFinite(metaMin) ? metaMin : Math.min(...meanScores);
  const dataMax = Number.isFinite(metaMax) ? metaMax : Math.max(...meanScores);
  const min = Math.floor(dataMin / CHRONOS_STEP) * CHRONOS_STEP;
  const max = Math.ceil(dataMax / CHRONOS_STEP) * CHRONOS_STEP;

  return {
    max: Math.max(max, min + CHRONOS_STEP),
    min,
  };
}

function syncScoreControls(event) {
  const min = Number(els.min.value);
  const max = Number(els.max.value);

  if (min > max && event.target === els.min) {
    els.max.value = min.toFixed(2);
  }

  if (max < min && event.target === els.max) {
    els.min.value = max.toFixed(2);
  }

  state.scoreMin = Number(els.min.value);
  state.scoreMax = Number(els.max.value);
  els.minValue.textContent = formatScore(state.scoreMin);
  els.maxValue.textContent = formatScore(state.scoreMax);
  render();
}

function syncChronosControls() {
  state.chronosMax = Number(els.chronosMax.value);
  renderChronosControls();
  render();
}

function renderChronosControls() {
  els.chronosMaxValue.textContent = formatChronosScore(state.chronosMax);
}

async function searchCurrentProtein() {
  const query = normalizeProtein(els.input.value || "TP53");
  els.input.value = query;
  state.center = query;
  state.selected = null;
  closeProteinCard();
  state.scoreMin = Number(els.min.value);
  state.scoreMax = Number(els.max.value);
  state.physicalOnly = els.physicalOnly.checked;

  setLoading(true);

  try {
    state.edges = await fetchStringInteractions(
      query,
      state.scoreMin,
      currentNetworkType(),
  );
    state.source = state.physicalOnly ? "STRING physical" : "STRING API";
  } catch (error) {
    console.warn("STRING API failed", error);
    state.edges = state.physicalOnly ? [] : findCachedInteractions(query);
    state.source = state.edges.length ? "local cache" : "no data";
  } finally {
    setLoading(false);
    render();
  }
}

async function fetchStringInteractions(protein, minScore, networkType) {
  const params = new URLSearchParams({
    identifier: protein,
    network_type: networkType,
    species: "9606",
    required_score: String(Math.round(minScore * 1000)),
    limit: String(STRING_LIMIT),
    caller_identity: "ppi-star-map-local",
  });
  const response = await fetch(`${STRING_API}?${params.toString()}`);

  if (!response.ok) {
    throw new Error(`STRING returned ${response.status}`);
  }

  return parseStringTsv(await response.text(), protein);
}

function currentNetworkType() {
  return state.physicalOnly ? "physical" : "functional";
}

function parseStringTsv(tsv, query) {
  const lines = tsv.trim().split(/\r?\n/).filter(Boolean);
  const headers = lines.shift()?.split("\t") || [];

  return lines
    .map((line) => {
      const row = Object.fromEntries(
        line.split("\t").map((value, index) => [headers[index], value]),
      );
      return {
        stringId_A: row.stringId_A,
        stringId_B: row.stringId_B,
        preferredName_A: row.preferredName_A,
        preferredName_B: row.preferredName_B,
        center: row.preferredName_A || query,
        partner: row.preferredName_B,
        partnerStringId: row.stringId_B,
        score: Number(row.score),
        evidence: Object.fromEntries(
          Object.keys(evidenceLabels).map((key) => [key, Number(row[key] || 0)]),
        ),
      };
    })
    .filter((edge) => edge.partner && Number.isFinite(edge.score))
    .sort((a, b) => b.score - a.score || a.partner.localeCompare(b.partner));
}

function findCachedInteractions(query) {
  return state.cache.interactions
    .filter(
      (edge) =>
        normalizeProtein(edge.preferredName_A) === query ||
        normalizeProtein(edge.preferredName_B) === query,
    )
    .map((edge) => {
      const centerIsA = normalizeProtein(edge.preferredName_A) === query;
      return {
        ...edge,
        center: centerIsA ? edge.preferredName_A : edge.preferredName_B,
        partner: centerIsA ? edge.preferredName_B : edge.preferredName_A,
        partnerStringId: centerIsA ? edge.stringId_B : edge.stringId_A,
      };
    })
    .sort((a, b) => b.score - a.score || a.partner.localeCompare(b.partner));
}

function render() {
  const edges = filteredEdges();
  const box = els.graphWrap.getBoundingClientRect();
  const viewportWidth = Math.max(520, box.width || 900);
  const viewportHeight = Math.max(520, box.height || 620);
  const layout = layoutConstellation(edges, viewportWidth, viewportHeight);

  els.graph.style.width = `${layout.width}px`;
  els.graph.style.height = `${layout.height}px`;
  els.graph.setAttribute("viewBox", `0 0 ${layout.width} ${layout.height}`);
  els.graph.replaceChildren();
  els.empty.hidden = edges.length > 0;
  els.title.textContent = `Potential effectors for ${state.center}`;
  els.visibleCount.textContent = String(edges.length);
  renderProteinCard(edges);

  if (!edges.length) {
    renderCenter(layout.center, state.center);
    centerGraph();
    return;
  }

  const edgeLayer = svg("g", { class: "edge-layer" });
  const nodeLayer = svg("g", { class: "node-layer" });
  const dustLayer = renderSpaceDust(layout);

  layout.nodes.forEach((node) => {
    edgeLayer.append(renderEdge(layout.center, node));
    nodeLayer.append(renderPartner(node));
  });

  els.graph.append(dustLayer, edgeLayer, nodeLayer);
  renderCenter(layout.center, state.center);
  centerGraph();
}

function filteredEdges() {
  const isChronosActive = chronosFilterActive();

  return state.edges.filter((edge) => {
    const scoreMatch = edge.score >= state.scoreMin && edge.score <= state.scoreMax;
    const chronosScore = getChronosSummary(edge.partner)?.mean_score;
    const chronosMatch =
      !isChronosActive ||
      (Number.isFinite(chronosScore) && chronosScore <= state.chronosMax);
    const essentialMatch =
      !state.onlyEssential || state.essential.has(normalizeProtein(edge.partner));
    const ligandMatch =
      !state.knownLigandOnly || hasKnownLigand(edge.partner);
    return scoreMatch && chronosMatch && essentialMatch && ligandMatch;
  });
}

function findVisibleEdge(edges, partner) {
  const symbol = normalizeProtein(partner);
  return edges.find((edge) => normalizeProtein(edge.partner) === symbol);
}

function hasKnownLigand(symbol) {
  return state.knownLigands.has(normalizeProtein(symbol));
}

function getChronosSummary(symbol) {
  return state.chronos.get(normalizeProtein(symbol));
}

function renderProteinCard(edges) {
  if (state.activeCard === "effector" && state.selected) {
    const edge = findVisibleEdge(edges, state.selected);

    if (edge) {
      renderEffectorCard(edge);
      return;
    }

    closeProteinCard();
  }

  renderTargetCard();
}

function renderTargetCard() {
  const symbol = normalizeProtein(state.center);
  const protein = state.ipdbProteins.get(symbol);
  const evidence = state.ipdbEvidence.get(symbol) || {};

  els.targetCard.classList.remove("is-effector-card");
  els.targetCard.setAttribute("aria-label", "Target protein card");
  els.targetName.textContent = protein?.name || symbol;
  els.targetSections.replaceChildren(
    renderCancerTypesSection(evidence.target_scores, evidence.expression),
    renderBindingSection(evidence.best_binding),
    renderLiteraturePocketsSection(evidence.literature_pockets),
    renderCompartmentsSection(evidence.compartments),
  );
}

function renderEffectorCard(edge) {
  const target = normalizeProtein(state.center);
  const effector = normalizeProtein(edge.partner);
  const protein = state.ipdbProteins.get(effector);
  const evidence = state.ipdbEvidence.get(effector) || {};

  els.targetCard.classList.add("is-effector-card");
  els.targetCard.setAttribute("aria-label", "Effector protein card");
  els.targetName.textContent = protein?.name || effector;
  els.targetSections.replaceChildren(
    renderStringSection(edge),
    renderDepMapSection(effector),
    renderSharedCompartmentsSection(effector, target),
    renderBindingSection(evidence.best_binding),
    renderLiteraturePocketsSection(evidence.literature_pockets),
    renderRiptacSection(edge),
  );
}

function renderStringSection(edge) {
  const section = cardSection(`STRING interaction with ${state.center}`);
  const channels = Object.entries(edge.evidence || {})
    .filter(([, value]) => Number(value) > 0)
    .sort(([, a], [, b]) => b - a)
    .slice(0, 2)
    .map(([key, value]) => `${evidenceLabels[key]} ${formatEvidence(value)}`);

  section.append(
    cardValue(`Confidence: ${edge.score.toFixed(3)}`),
    cardMeta(channels.join(" \u00b7 ")),
  );
  return section;
}

function renderDepMapSection(symbol) {
  const section = cardSection("DepMap essentiality");
  const summary = getChronosSummary(symbol);

  if (!summary) {
    section.append(cardValue("No DepMap evidence"));
    return section;
  }

  section.append(
    cardValue(`Mean Chronos: ${formatChronosScore(summary.mean_score)}`),
    cardMeta(
      `${formatPercent(summary.dependency_fraction)} cell lines <= ${formatChronosScore(
        state.dependencyScoreThreshold,
      )}`,
    ),
  );
  return section;
}

function renderSharedCompartmentsSection(effector, target) {
  const section = cardSection("Shared compartments");
  const effectorCompartments = knownCompartments(
    state.ipdbEvidence.get(effector)?.compartments,
  );
  const targetCompartments = knownCompartments(
    state.ipdbEvidence.get(target)?.compartments,
  );

  if (!effectorCompartments.length || !targetCompartments.length) {
    section.append(cardValue("Insufficient compartment annotation"));
    return section;
  }

  const targetSet = new Set(targetCompartments);
  const shared = effectorCompartments.filter((compartment) =>
    targetSet.has(compartment),
  );

  if (!shared.length) {
    section.append(cardValue("No shared compartment"));
    return section;
  }

  section.append(renderCompartmentList(shared));
  return section;
}

function renderRiptacSection(edge) {
  const section = cardSection("Target-first RIPTAC");
  const bestRiptac = getTargetFirstRiptacRows(edge)[0];

  section.classList.add("is-wide");

  if (!bestRiptac) {
    section.append(
      cardValue(`No ${edge.partner} -> ${state.center} annotation`),
    );
    return section;
  }

  section.append(
    cardValue(`IPDB pair score: ${bestRiptac.score.toFixed(3)}`),
    cardMeta(
      `${bestRiptac.left} -> ${bestRiptac.right} \u00b7 ${formatCancerType(
        bestRiptac.cancer_type,
      )}`,
    ),
  );
  return section;
}

function renderCancerTypesSection(targetScores = [], expression) {
  const section = cardSection("Top 3 cancer types");
  const scores = Array.isArray(targetScores)
    ? targetScores
        .filter((targetScore) => Number.isFinite(Number(targetScore.score)))
        .slice(0, 3)
    : [];

  section.classList.add("is-wide");

  if (scores.length) {
    section.append(renderTargetScoreList(scores));
    return section;
  }

  if (expression) {
    section.append(
      cardValue(formatCancerType(expression.cancer_type)),
      cardMeta(formatExpressionMeta(expression)),
    );
    return section;
  }

  section.append(cardValue("No IPDB cancer-type evidence"));
  return section;
}

function renderTargetScoreList(scores) {
  const list = document.createElement("ol");
  list.className = "card-metric-list";
  list.replaceChildren(
    ...scores.map((targetScore) => {
      const item = document.createElement("li");
      const cancerType = document.createElement("span");
      const score = document.createElement("strong");

      cancerType.textContent = formatCancerType(targetScore.cancer_type);
      score.textContent = formatMetric(targetScore.score);
      item.append(cancerType, score);

      return item;
    }),
  );

  return list;
}

function renderLiteraturePocketsSection(literaturePockets = []) {
  const section = cardSection("Literature pockets");
  const pockets = Array.isArray(literaturePockets)
    ? literaturePockets.slice(0, 3)
    : [];

  if (!pockets.length) {
    section.append(cardValue("No literature pocket evidence"));
    return section;
  }

  const list = document.createElement("ol");
  list.className = "card-metric-list";
  list.replaceChildren(
    ...pockets.map((pocket) => {
      const item = document.createElement("li");
      const label = document.createElement("span");
      const ligandabilityScore = Number(pocket.ligandability_score);

      label.textContent =
        pocket.evidence_label || formatLabel(pocket.pocket_type);
      item.append(label);

      if (Number.isFinite(ligandabilityScore)) {
        const score = document.createElement("strong");
        score.textContent = formatMetric(ligandabilityScore);
        item.append(score);
      }

      return item;
    }),
  );

  section.append(list);
  return section;
}

function formatExpressionMeta(expression) {
  const score = Number(expression.overexpression_score);
  const foldChange = Number(expression.log2_fold_change);
  return [
    Number.isFinite(score) ? `Overexpression score: ${formatMetric(score)}` : null,
    expression.best_assay ? String(expression.best_assay).toUpperCase() : null,
    Number.isFinite(foldChange) ? `log2FC ${formatMetric(foldChange)}` : null,
  ].filter(Boolean).join(" \u00b7 ");
}

function renderBindingSection(binding) {
  const section = cardSection("Best BindingDB binder");

  if (!binding) {
    section.append(cardValue("No BindingDB hit"));
    return section;
  }

  const value = Number(binding.value_nm);
  const affinity = Number.isFinite(value) ? `${formatNm(value)} nM` : null;
  const meta = [affinity, binding.assay_type].filter(Boolean);

  section.append(
    cardValue(binding.ligand_name || binding.ligand_id || "BindingDB hit"),
    cardMeta(meta.join(" \u00b7 ")),
  );
  return section;
}

function renderCompartmentsSection(compartments = []) {
  const known = knownCompartments(compartments);
  const section = cardSection("Compartments");

  if (!known.length) {
    section.append(cardValue("No compartment annotation"));
    return section;
  }

  section.append(renderCompartmentList(known));
  return section;
}

function renderCompartmentList(compartments) {
  const chips = document.createElement("div");
  chips.className = "compartment-list";
  chips.replaceChildren(
    ...compartments.map((compartment) => {
      const chip = document.createElement("span");
      chip.textContent = formatLabel(compartment);
      return chip;
    }),
  );

  return chips;
}

function cardSection(label) {
  const section = document.createElement("div");
  section.className = "target-card-section";
  section.append(cardLabel(label));
  return section;
}

function cardLabel(text) {
  const label = document.createElement("span");
  label.className = "target-card-label";
  label.textContent = text;
  return label;
}

function cardValue(text) {
  const value = document.createElement("strong");
  value.textContent = text;
  return value;
}

function cardMeta(text) {
  const meta = document.createElement("small");
  meta.textContent = text;
  meta.hidden = !text;
  return meta;
}

function chronosFilterActive() {
  return state.chronosMax < state.chronosBounds.max - 0.0001;
}

function renderSpaceDust(layout) {
  const layer = svg("g", {
    class: "space-dust",
    "aria-hidden": "true",
  });
  const area = layout.width * layout.height;
  const starCount = Math.round(clamp(area / 6200, 90, 360));

  layer.append(
    svg("ellipse", {
      class: "nebula-cloud cloud-a",
      cx: layout.center.x - layout.width * 0.08,
      cy: layout.center.y + layout.height * 0.02,
      rx: layout.width * 0.34,
      ry: layout.height * 0.24,
    }),
    svg("ellipse", {
      class: "nebula-cloud cloud-b",
      cx: layout.center.x + layout.width * 0.14,
      cy: layout.center.y - layout.height * 0.1,
      rx: layout.width * 0.22,
      ry: layout.height * 0.16,
    }),
  );

  for (let index = 0; index < starCount; index += 1) {
    const seed = `${state.center}:${layout.width}:${layout.height}:${index}`;
    const twinkle = stableNoise(`${seed}:twinkle`);
    layer.append(
      svg("circle", {
        class: `dust-star ${twinkle > 0.92 ? "bright" : ""}`,
        cx: stableNoise(`${seed}:x`) * layout.width,
        cy: stableNoise(`${seed}:y`) * layout.height,
        r: 0.5 + stableNoise(`${seed}:r`) * 1.35,
        style: `--twinkle:${0.24 + twinkle * 0.68}`,
      }),
    );
  }

  return layer;
}

function renderEdge(center, node) {
  const dx = node.x - center.x;
  const dy = node.y - center.y;
  const bend = (stableNoise(`${node.edge.partner}:bend`) - 0.5) * node.radius * 0.16;
  const midX = center.x + dx * 0.54 - Math.sin(node.angle) * bend;
  const midY = center.y + dy * 0.54 + Math.cos(node.angle) * bend;

  return svg("path", {
    class: `edge ${node.tier}`,
    d: `M ${center.x.toFixed(2)} ${center.y.toFixed(2)} Q ${midX.toFixed(
      2,
    )} ${midY.toFixed(2)} ${node.x.toFixed(2)} ${node.y.toFixed(2)}`,
    style: `--edge-opacity:${(0.13 + node.strength * 0.55).toFixed(
      3,
    )}; --edge-width:${(0.7 + node.strength * 3.4).toFixed(2)}px;`,
  });
}

function renderPartner(node) {
  const isEssential = state.essential.has(normalizeProtein(node.edge.partner));
  const riptacRows = getTargetFirstRiptacRows(node.edge);
  const bestRiptac = riptacRows[0];
  const chronosSummary = getChronosSummary(node.edge.partner);
  const chronosLabel = chronosSummary
    ? `, mean Chronos ${formatChronosScore(chronosSummary.mean_score)}`
    : "";
  const riptacLabel = bestRiptac
    ? `, RIPTAC ${bestRiptac.score.toFixed(3)} for ${bestRiptac.left} to ${
        bestRiptac.right
      }`
    : "";
  const group = svg("g", {
    class: `node ${node.tier}${node.hasLabel ? " has-label" : " is-particle"}${
      isEssential ? " is-essential" : ""
    }${bestRiptac ? " has-riptac" : ""}${
      state.selected === node.edge.partner ? " selected" : ""
    }`,
    transform: `translate(${node.x} ${node.y})`,
    "data-partner": node.edge.partner,
    tabindex: "0",
    role: "button",
    "aria-expanded": String(
      state.activeCard === "effector" && state.selected === node.edge.partner,
    ),
    "aria-label": `${node.edge.partner}, score ${node.edge.score.toFixed(
      3,
    )}${chronosLabel}${riptacLabel}`,
  });

  group.append(
    svg("title", {
      text: `${node.edge.partner} - STRING ${node.edge.score.toFixed(
        3,
      )}${chronosLabel.replace(",", " -")}${riptacLabel.replace(",", " -")}`,
    }),
  );

  if (node.hasLabel) {
    group.append(
      svg("ellipse", {
        class: "node-body",
        rx: node.rx,
        ry: node.ry,
        fill: node.fill,
        stroke: node.stroke,
        "stroke-width": node.strokeWidth,
        style: node.style,
      }),
      svg("text", {
        "dominant-baseline": "central",
        "font-size": node.font.toFixed(1),
        text: node.edge.partner,
      }),
    );
  } else {
    group.append(
      svg("circle", {
        class: "node-dot",
        r: node.rx,
        fill: node.fill,
        stroke: node.stroke,
        "stroke-width": node.strokeWidth,
        style: node.style,
      }),
      svg("g", {
        class: "node-callout",
        transform: `translate(0 ${-node.calloutOffset})`,
      }),
    );

    group.querySelector(".node-callout").append(
      svg("rect", {
        x: -node.calloutWidth / 2,
        y: -node.calloutHeight / 2,
        width: node.calloutWidth,
        height: node.calloutHeight,
        rx: 6,
      }),
      svg("text", {
        "dominant-baseline": "central",
        "font-size": node.calloutFont.toFixed(1),
        text: node.edge.partner,
      }),
    );
  }

  if (isEssential) {
    group.append(
      svg("ellipse", {
        class: "essential-ring",
        rx: node.rx + 6,
        ry: node.ry + 6,
      }),
    );
  }

  if (bestRiptac) {
    group.append(
      svg("ellipse", {
        class: "riptac-ring",
        rx: node.rx + 10,
        ry: node.ry + 10,
      }),
    );
  }

  return group;
}

function renderCenter(center, label) {
  const group = svg("g", {
    class: "center-node",
    transform: `translate(${center.x} ${center.y})`,
    tabindex: "0",
    role: "button",
    "aria-expanded": String(state.activeCard === "target"),
    "aria-label": `Open target protein card for ${label}`,
  });
  group.append(
    svg("title", { text: `Open target protein card for ${label}` }),
    svg("circle", { class: "center-orbit", r: 76 }),
    svg("circle", { class: "center-halo", r: 56 }),
    svg("circle", { class: "center-core", r: 42 }),
    svg("text", {
      "dominant-baseline": "central",
      "font-size": 18,
      text: label,
    }),
  );
  els.graph.append(group);
}

function getTargetFirstRiptacRows(edge) {
  const effector = normalizeProtein(edge.partner);
  const target = normalizeProtein(edge.center || state.center);
  const pair = state.ipdbPairs.get(normalizePairKey(effector, target));

  return (pair?.riptac || [])
    .filter(
      (row) =>
        normalizeProtein(row.left) === effector &&
        normalizeProtein(row.right) === target,
    )
    .sort((a, b) => b.score - a.score);
}

function layoutConstellation(edges, viewportWidth, viewportHeight) {
  const center = { x: viewportWidth / 2, y: viewportHeight / 2 };

  if (!edges.length) {
    return { center, height: viewportHeight, nodes: [], width: viewportWidth };
  }

  const scores = edges.map((edge) => edge.score);
  const minScore = Math.min(...scores);
  const maxScore = Math.max(...scores);
  const crowding = clamp((edges.length - 56) / 320, 0, 1);
  const labelBudget =
    edges.length <= 56 ? edges.length : Math.round(42 - crowding * 20);
  const shells = [];
  let cursor = 0;
  let radius = 122;
  const nodes = edges.map((edge, index) =>
    createNode(
      edge,
      visualStrength(edge.score, index, edges.length, minScore, maxScore),
      index,
      labelBudget,
      crowding,
    ),
  );

  while (cursor < nodes.length) {
    const shell = [];
    const circumference = Math.PI * 2 * radius;
    const usableCircumference = circumference * (0.86 + crowding * 0.08);
    let usedCircumference = 0;

    while (cursor < nodes.length) {
      const node = nodes[cursor];
      const slotWidth = node.rx * 1.85 + 18;

      if (
        shell.length >= 3 &&
        usedCircumference + slotWidth > usableCircumference
      ) {
        break;
      }

      shell.push(node);
      usedCircumference += slotWidth;
      cursor += 1;
    }

    shells.push({ nodes: shell, radius });
    radius += Math.max(
      34,
      22 + Math.max(...shell.map((node) => node.ry * 2.2)),
    );
  }

  const placedNodes = [];

  shells.forEach((shell, shellIndex) => {
    const shellJitter = Math.min(
      12,
      (shells[shellIndex + 1]?.radius - shell.radius || 64) * 0.16,
    );
    const startAngle =
      -Math.PI / 2 +
      shellIndex * 0.37 +
      (stableNoise(`${state.center}:shell:${shellIndex}`) - 0.5) * 0.54;
    const direction = shellIndex % 2 === 0 ? 1 : -1;

    shell.nodes.forEach((node, index) => {
      const angle =
        startAngle +
        direction * (index / shell.nodes.length) * Math.PI * 2 +
        (stableNoise(`${node.edge.partner}:angle`) - 0.5) * 0.08;
      const radialJitter = (stableNoise(node.edge.partner) - 0.5) * shellJitter;

      node.angle = angle;
      node.radius = shell.radius + radialJitter;
      node.x = Math.cos(angle) * node.radius;
      node.y = Math.sin(angle) * node.radius;
      placedNodes.push(node);
    });
  });

  const outerRadius =
    Math.max(
      ...placedNodes.map(
        (node) => Math.hypot(node.x, node.y) + Math.max(node.rx, node.ry),
      ),
    ) + 84;
  const diameter = Math.ceil(outerRadius * 2);
  const width = Math.max(viewportWidth, diameter);
  const height = Math.max(viewportHeight, diameter);
  const graphCenter = { x: width / 2, y: height / 2 };

  placedNodes.forEach((node) => {
    node.x += graphCenter.x;
    node.y += graphCenter.y;
  });

  return { center: graphCenter, height, nodes: placedNodes, width };
}

function createNode(edge, strength, index, labelBudget, crowding) {
  const hasLabel = index < labelBudget;
  const tier = strength > 0.72 ? "core" : strength > 0.38 ? "mid" : "dust";
  const nodeScale = 1 - crowding * 0.34;
  const cyan = Math.round(146 + strength * 82);
  const font = (9.5 + strength * 7.7) * Math.max(0.82, nodeScale);
  const calloutFont = 11 + strength * 2;
  const labelWidth = estimateLabelWidth(edge.partner, font) + 30;
  const calloutWidth = estimateLabelWidth(edge.partner, calloutFont) + 28;
  const dotRadius = (4.2 + strength * 6.8) * Math.max(0.78, nodeScale);

  return {
    calloutFont,
    calloutHeight: 28,
    calloutOffset: 24 + dotRadius,
    calloutWidth: Math.max(58, calloutWidth),
    edge,
    fill: `rgba(${Math.round(18 + strength * 52)}, ${cyan}, 255, ${
      0.35 + strength * 0.42
    })`,
    font,
    hasLabel,
    rx: hasLabel ? Math.max(28, (labelWidth * nodeScale) / 2) : dotRadius,
    ry: hasLabel ? (15 + strength * 9.5) * nodeScale : dotRadius,
    strength,
    stroke: `rgba(224, 252, 255, ${0.4 + strength * 0.56})`,
    strokeWidth: hasLabel ? 0.9 + strength * 1.4 : 0.8 + strength * 0.9,
    style: `--node-glow: rgba(142, 246, 255, ${
      0.16 + strength * 0.62
    }); --node-glow-size: ${(5 + strength * 21).toFixed(1)}px;`,
    tier,
    x: 0,
    y: 0,
  };
}

function estimateLabelWidth(label, fontSize) {
  return label.length * fontSize * 0.64;
}

function visualStrength(score, index, total, minScore, maxScore) {
  const scoreStrength =
    maxScore === minScore ? 1 : (score - minScore) / (maxScore - minScore);
  const rankStrength = total <= 1 ? 1 : 1 - index / (total - 1);

  return Math.max(0.12, Math.min(1, scoreStrength * 0.78 + rankStrength * 0.22));
}

function selectNode(event) {
  if (event.type === "keydown" && !["Enter", " "].includes(event.key)) return;

  const centerNode = event.target.closest(".center-node");
  if (centerNode) {
    event.preventDefault();
    clearSelectedNodes();
    openProteinCard("target");
    renderTargetCard();
    return;
  }

  const node = event.target.closest(".node");
  if (!node) return;

  event.preventDefault();
  const edge = findVisibleEdge(filteredEdges(), node.dataset.partner);
  if (!edge) return;

  clearSelectedNodes();
  state.selected = node.dataset.partner;
  openProteinCard("effector");
  renderEffectorCard(edge);
  node.classList.add("selected");
}

function openProteinCard(type) {
  state.activeCard = type;
  els.targetCard.hidden = false;
  syncCardExpansion();
}

function closeProteinCard() {
  clearSelectedNodes();
  state.activeCard = null;
  els.targetCard.hidden = true;
  syncCardExpansion();
}

function clearSelectedNodes() {
  state.selected = null;
  els.graph
    .querySelectorAll(".node.selected")
    .forEach((selectedNode) => selectedNode.classList.remove("selected"));
}

function syncCardExpansion() {
  els.graph
    .querySelector(".center-node")
    ?.setAttribute("aria-expanded", String(state.activeCard === "target"));
  els.graph.querySelectorAll(".node").forEach((node) => {
    node.setAttribute(
      "aria-expanded",
      String(
        state.activeCard === "effector" &&
          state.selected === node.dataset.partner,
      ),
    );
  });
}

function centerGraph() {
  requestAnimationFrame(() => {
    els.graphWrap.scrollLeft =
      (els.graphWrap.scrollWidth - els.graphWrap.clientWidth) / 2;
    els.graphWrap.scrollTop =
      (els.graphWrap.scrollHeight - els.graphWrap.clientHeight) / 2;
  });
}

function stableNoise(value) {
  let hash = 2166136261;

  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 16777619);
  }

  return (hash >>> 0) / 4294967295;
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function setLoading(isLoading) {
  els.loading.hidden = !isLoading;
  els.apply.disabled = isLoading;
  els.form.querySelector("button").disabled = isLoading;
}

function normalizeProtein(value) {
  return String(value || "")
    .trim()
    .toUpperCase();
}

function normalizeEntryMap(entries) {
  return new Map(
    Object.entries(entries).map(([symbol, value]) => [
      normalizeProtein(symbol),
      value,
    ]),
  );
}

function knownCompartments(compartments = []) {
  return compartments.filter(
    (compartment) => compartment && compartment !== "unknown",
  );
}

function normalizePairKey(left, right) {
  return [normalizeProtein(left), normalizeProtein(right)].sort().join("::");
}

function formatScore(value) {
  return value.toFixed(2);
}

function formatChronosScore(value) {
  return Number(value).toFixed(2);
}

function formatEvidence(value) {
  return Number(value).toFixed(2);
}

function formatMetric(value) {
  return Number(value).toFixed(2);
}

function formatNm(value) {
  if (value < 0.01) return value.toPrecision(1);
  if (value < 1) return value.toPrecision(2);
  if (value < 10) return value.toFixed(2);
  if (value < 100) return value.toFixed(1);
  return value.toLocaleString(undefined, { maximumFractionDigits: 0 });
}

function formatPercent(value) {
  return `${(Number(value) * 100).toFixed(1)}%`;
}

function formatCancerType(value) {
  return formatLabel(value || "Unknown cancer type");
}

function formatLabel(value) {
  return String(value || "").replaceAll("_", " ");
}

function svg(tag, attrs = {}) {
  const el = document.createElementNS(SVG_NS, tag);

  Object.entries(attrs).forEach(([key, value]) => {
    if (key === "text") {
      el.textContent = value;
    } else {
      el.setAttribute(key, value);
    }
  });

  return el;
}
