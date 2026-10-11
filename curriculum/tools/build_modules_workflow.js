export const meta = {
  name: 'build-modules',
  description: 'Write or complete curriculum modules from briefs, adversarially verify the math, fix findings',
  phases: [
    { title: 'Write', detail: 'one writer per module: complete the four files from the brief, script must pass' },
    { title: 'Verify', detail: 'one adversarial math reviewer per module: derivations, numbers, consistency, style' },
    { title: 'Fix', detail: 'apply verified findings and re-run the script' },
  ],
}

const ROOT = '/home/user/sinjinbhatia-hub/curriculum'
const modules = (args && args.modules) || []
if (!modules.length) throw new Error('args.modules required: [{code, dir, brief, mode?}]')

const WRITE_SCHEMA = {
  type: 'object',
  properties: {
    status: { type: 'string', enum: ['complete', 'partial', 'failed'] },
    script_passes: { type: 'boolean' },
    runtime_s: { type: 'number' },
    word_counts: { type: 'object', properties: { notes: { type: 'number' }, exercises: { type: 'number' }, solutions: { type: 'number' } } },
    deviations: { type: 'string' },
  },
  required: ['status', 'script_passes'],
}

const VERIFY_SCHEMA = {
  type: 'object',
  properties: {
    script_passes: { type: 'boolean' },
    issues: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          file: { type: 'string' },
          location: { type: 'string' },
          description: { type: 'string' },
          suggested_fix: { type: 'string' },
        },
        required: ['severity', 'file', 'location', 'description', 'suggested_fix'],
      },
    },
    summary: { type: 'string' },
  },
  required: ['script_passes', 'issues', 'summary'],
}

const FIX_SCHEMA = {
  type: 'object',
  properties: {
    fixed: { type: 'number' },
    skipped: { type: 'number' },
    script_passes: { type: 'boolean' },
    notes: { type: 'string' },
  },
  required: ['fixed', 'skipped', 'script_passes'],
}

function writePrompt(m) {
  return `You are completing one module of a self-study curriculum.

1. Read ${ROOT}/research/WRITER-BRIEF.md and follow it exactly (it tells you to read STYLE.md, TEMPLATE.md, and the exemplar module).
2. Read your module brief at ${ROOT}/research/briefs/${m.brief}. It states the module code, folder, topics, required derivations, worked examples, project contents, and (if present) a STATUS/REPAIR section describing which files already exist and what is missing or broken; the folder's actual contents take precedence over that section if they differ.
3. Inspect the module folder ${ROOT}/${m.dir}/ first. Files that already exist were written by a previous writer that was interrupted: read them in full, keep what is good, and complete or repair rather than starting over, unless a file is clearly truncated garbage. Never write outside this folder.
4. Produce the four files (notes.md, exercises.md, solutions.md, project.py) satisfying the brief and the template. Run \`cd ${ROOT}/${m.dir} && MPLBACKEND=Agg python3 project.py\`; it must exit 0 in under ~90 s using only numpy/scipy/matplotlib/sympy (scikit-learn or torch only where the brief says so). Fix until it passes. Every number quoted in a worked example or solution must be asserted or printed by the script.
5. Self-check against STYLE.md's "things that get a module sent back" list.
6. Return the structured result. In \`deviations\`, say in under 120 words where you departed from the brief and why. Do not paste file contents.`
}

function verifyPrompt(m) {
  return `You are an adversarial mathematical reviewer for one module of a self-study curriculum. Your job is to find errors, not to praise. Assume there are mistakes and hunt for them.

Module folder: ${ROOT}/${m.dir}/ (files notes.md, exercises.md, solutions.md, project.py). Brief: ${ROOT}/research/briefs/${m.brief} (if present). Style rules: ${ROOT}/STYLE.md. Template: ${ROOT}/TEMPLATE.md.

Do all of the following and report every problem you find:
A. MATH. Read notes.md section 3 and section 4 and solutions.md line by line. For every derivation, check each numbered step follows from the previous one; check signs, factors of 2 and pi, transposes vs inverses, conjugates, index ranges, and hypotheses actually used. For every theorem statement, check the hypotheses are sufficient and the proof uses them. For every numerical value (worked examples and solutions), recompute it yourself with python3 (numpy/sympy/scipy) in a scratch script under /tmp/claude-0 (not in the module folder) and compare. Flag any value that disagrees beyond rounding.
B. CONSISTENCY. Check that exercises.md and solutions.md have the same numbering and that every exercise has a complete solution with no "similarly", "left to the reader", or "it can be shown" for anything central. Check that the conventions box at the top of notes.md matches what the notes actually do (e.g. which slot of a complex inner product is conjugated; where 2 pi sits in the Fourier convention; qubit ordering). Check that cross-references to other modules (e.g. "LA-8 Theorem 3.4") name an object that exists in that module's notes if that module exists under ${ROOT}; if the referenced module does not exist yet, only flag it if the claim about it is wrong in substance.
C. CODE. Run \`cd ${ROOT}/${m.dir} && MPLBACKEND=Agg timeout 300 python3 project.py\`. Report whether it exits 0 and the runtime. Check that the asserts correspond to the claims in the notes (an assert that tests something weaker than the notes claim is an issue). Check for imports outside numpy/scipy/matplotlib/sympy (scikit-learn/torch only if the brief allows).
D. TEMPLATE/STYLE. Sections 0 to 10 present with the template headings; conventions box; "Try it first"/"What just happened" around key derivations; at least three worked examples; A/B/C exercises; section 8 is specific, not generic; section 10 agrees with the brief's pairing list.
E. PEDAGOGY (lower severity). Places where section 2's intuition is just section 3 with fewer symbols; a derivation that skips the step a learner would get stuck on; invented or dubious historical attributions (flag any you cannot confirm from your own knowledge as medium severity with the wording "unverified attribution").

Severity: high = a mathematical error, a wrong number, a failing or weakened script, a missing file or missing required section, a contradiction between files; medium = an incomplete derivation or solution, a convention inconsistency, an unverified attribution; low = pedagogy and wording.
Be specific: file, exact location (section/theorem/exercise number or a short quoted phrase), what is wrong, and the concrete fix. Do not edit any files yourself. Return the structured result; put the count of issues by severity in \`summary\`.`
}

function fixPrompt(m, v) {
  return `You are fixing verified review findings in one module of a self-study curriculum.

Module folder: ${ROOT}/${m.dir}/. Style rules: ${ROOT}/STYLE.md. Brief: ${ROOT}/research/briefs/${m.brief} (if present).

The reviewer's findings (JSON):
${JSON.stringify(v.issues, null, 1)}

Reviewer summary: ${v.summary}

Instructions:
1. Fix every high- and medium-severity finding. For each, first confirm it is real (re-derive or recompute); if a finding is wrong, skip it and say why in \`notes\`. Fix low-severity findings when the fix is local and clearly an improvement.
2. Keep notes.md, exercises.md, solutions.md, and project.py mutually consistent: if you change a number or a claim in one file, change it everywhere it appears, and make sure project.py asserts the corrected value.
3. Re-run \`cd ${ROOT}/${m.dir} && MPLBACKEND=Agg timeout 300 python3 project.py\` until it exits 0.
4. Re-read each passage you edited once more for correctness. Do not write outside the module folder. Do not paste file contents in your result.
Return the structured result with counts of fixed and skipped findings and a short \`notes\` field (under 150 words).`
}

const results = await pipeline(
  modules,
  m => (m.mode === 'verify')
    ? Promise.resolve({ status: 'complete', script_passes: true, skipped_write: true })
    : agent(writePrompt(m), { label: `write:${m.code}`, phase: 'Write', schema: WRITE_SCHEMA }),
  (w, m) => {
    if (!w) { log(`${m.code}: writer died; skipping verify`); return null }
    return agent(verifyPrompt(m), { label: `verify:${m.code}`, phase: 'Verify', schema: VERIFY_SCHEMA }).then(v => ({ w, v }))
  },
  (r, m) => {
    if (!r || !r.v) { log(`${m.code}: no verification result`); return r }
    const actionable = r.v.issues.filter(i => i.severity !== 'low')
    if (!actionable.length && r.v.script_passes) { log(`${m.code}: verified clean (${r.v.issues.length} low-severity notes)`); return { ...r, f: { fixed: 0, skipped: 0, script_passes: true } } }
    return agent(fixPrompt(m, r.v), { label: `fix:${m.code}`, phase: 'Fix', schema: FIX_SCHEMA }).then(f => ({ ...r, f }))
  },
)

const summary = modules.map((m, i) => {
  const r = results[i]
  if (!r) return { code: m.code, status: 'FAILED (writer or pipeline died)' }
  return {
    code: m.code,
    write: r.w ? (r.w.skipped_write ? 'pre-existing' : `${r.w.status}, script ${r.w.script_passes ? 'passes' : 'FAILS'}`) : 'n/a',
    verify: r.v ? `${r.v.issues.length} issues (${r.v.issues.filter(i => i.severity === 'high').length} high, ${r.v.issues.filter(i => i.severity === 'medium').length} medium); script ${r.v.script_passes ? 'passes' : 'FAILS'}` : 'no result',
    fix: r.f ? `${r.f.fixed} fixed, ${r.f.skipped} skipped; script ${r.f.script_passes ? 'passes' : 'FAILS'}${r.f.notes ? '; ' + r.f.notes : ''}` : 'no result',
    deviations: r.w && r.w.deviations ? r.w.deviations : '',
  }
})
return summary
