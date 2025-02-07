/*
data attribute creates object that is bound to this

computed - computed property but with caching
method - computed property - invoked each time

v-model - links inputs to vue js data


https://codepen.io/pespantelis/pen/ojwgPB
https://www.raymondcamden.com/2018/02/08/building-table-sorting-and-pagination-in-vuejs

*/

/* EventBus is used to pass changes between components */
// const EventBus = new Vue();

Vue.filter("decimalPlaces", (value, num = 2) => {
    if (value == null) {
        return "N/A";
    } else {
        return value.toFixed(num);
    }
});

Vue.component('hklinfo-table', {
    data: function () {
        return {
            hklinfo: this.$root.hklinfo,
        }
    },
    mounted() {
        let show = this.$root.hklinfo ? true : false;
        toggleDisplay("hklinfo-title", show);
        window.collapseDiv('collapsible1');

      },
    template: `
    <div id="div1" class="content">
        <div id="table-container"> 
            <table id="hkl_table">
                <thead>
                        <th title="Name of, and link to, the crystallographic data file">Name</th>
                        <th title="Highest resolution of the crystallographic data">Resolution</th>
                        <th title="The space group of the crystallographic data">Space Group</th>
                        <th title="Inidicates the pressence of Non-crystallographic symmetry (calculated by CTRUNCATE)">Has NCS?</th>
                        <th title="Inidicates the pressence of twinning (calculated by CTRUNCATE)">Has Twinning?</th>
                        <th title="Inidicates the pressence of anisotropy (calculated by CTRUNCATE)">Has Anisotropy?</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>
                            <a v-if="hklinfo.hklin" :href="hklinfo.hklin">{{ hklinfo.name }}</a>
                            <a v-else>{{ hklinfo.name }}</a>
                        </td>
                        <td>{{ hklinfo.resolution | decimalPlaces }}</td>
                        <td>{{ hklinfo.space_group }}</td>
                        <td>{{ hklinfo.has_ncs }}</td>
                        <td>{{ hklinfo.has_twinning }}</td>
                        <td>{{ hklinfo.has_anisotropy }}</td>
                    </tr>
                </tbody>
            </table>
        </div>
    </div>`
});

Vue.component('pdb-table', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            homologs: this.$root.homologs,
        }
    },
    mounted() {
        let show = Object.keys(this.$root.homologs).length !== 0;
        toggleDisplay("pdb-title", show);
        toggleDisplay("div3", show);
      },
    template: `
    <div id="pdb-table-container"> 
        <table id="pdb_table">
            <thead>
                <tr>
                    <th title="Name of the homologue (&lt;PDB&gt;_&lt;CHAIN_ID&gt;)" 
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 0)">Name ↓</th>
                    <th title="PDB code of the homlogue" 
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 1)">PDB ↓</th>
                    <th title="Resolution of the Xtal data for the PDB entry" 
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 2)">Resolution ↓</th>
                    <th title="The region the homologue covers in the sequence"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 3)">Region ↓</th>
                    <th title="Start-stop of the homologue in the sequence"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 4)">Range ↓</th>
                    <th title="Length of the homologue"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 5)">Length ↓</th>
                    <th title="The estimated Log Liklihood Gain when the homologue is used in molecular replacement"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 6)">eLLG ↓</th>
                    <th title="The molecular weight of the homologue"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 7)">Mol. Wt. ↓</th>
                    <th title="The estimated Root Mean Square Deviation of the homologue to the target protein"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 8)">eRMSD ↓</th>
                    <th title="The sequence identity of the homologue to the target protein"
                        onclick="sortTableAndFeatures('pdb_table', sequence, ft2, PDB_rowFeatureMap, '#div3', 9)">Seq. Ident. ↓</th>
                </tr>
            </thead>
            <tbody>
                <tr v-for="homolog in homologs" :data-feature-name="'PDB-' + homolog.name">
                    <td>
                    <a v-if="homolog.pdb_file" :href="homolog.pdb_file">{{ homolog.name }}</a>
                    <a v-else>{{ homolog.name }}</a>
                    </td>
                    <td><a v-bind:href="homolog.pdb_url" target="_blank">{{ homolog.pdb_id }}</a></td>
                    <td>{{ homolog.resolution  | decimalPlaces }}</td>
                    <td>{{ homolog.region_id }}</td>
                    <td>{{ homolog.range }}</td>
                    <td>{{ homolog.length }}</td>
                    <td>{{ homolog.ellg }}</td>
                    <td>{{ homolog.molecular_weight }}</td>
                    <td>{{ homolog.rmsd }}</td>
                    <td>{{ homolog.seq_ident }}</td>
                </tr>
            </tbody>
        </table>
    </div>
    `
});

Vue.component('af-table', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            af_models: this.$root.af_models,
        }
    },
    mounted() {
        let show = Object.keys(this.$root.af_models).length !== 0;
        toggleDisplay("af-title", show);
        toggleDisplay("div4", show);
      },
    template: `
    <div id="af-table-container"> 
        <table id="af_table">
            <thead>
                <tr>
                    <th title="The name of the model in the AlphaFold database"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 0)">Name ↓</th>
                    <th title="The date the model was made"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 1)">Date made ↓</th>
                    <th title="Start-stop of the model in the sequence"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 2)">Range ↓</th>
                    <th title="Length of the model"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 3)">Length ↓</th>
                    <th title="The average pLDDT score of the model"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 4)">Average pLDDT ↓</th>
                    <th title="The H-score of the model"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 5)">H-score ↓</th>
                    <th title="The sequence identity of the model to the target protein"
                        onclick="sortTableAndFeatures('af_table', sequence, ft3, AF_rowFeatureMap, '#div4', 6)">Sequence Identity ↓</th>
                </tr>
            </thead>
            <tbody>
                <tr v-for="afmodel in af_models" :data-feature-name="afmodel.model_id">
                    <td><a v-bind:href="afmodel.model_url" target="_blank">{{ afmodel.model_id }}</a></td>
                    <td>{{ afmodel.date_made }}</td>
                    <td>{{ afmodel.range }}</td>
                    <td>{{ afmodel.length }}</td>
                    <td>{{ afmodel.avg_plddt | decimalPlaces }}</td>
                    <td>{{ afmodel.h_score }}</td>
                    <td>{{ afmodel.seq_ident }}</td>
                </tr>
            </tbody>
        </table>
    </div>
    `
});

Vue.component('bfvd-table', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            bfvd_models: this.$root.bfvd_models,
        }
    },
    mounted() {
        let show = Object.keys(this.$root.bfvd_models).length !== 0;
        toggleDisplay("bfvd-title", show);
        toggleDisplay("div5", show);
      },
    template: `
    <div id="bfvd-table-container"> 
        <table id="bfvd_table">
            <thead>
                <tr>
                    <th title="The name of the model in the Big Fantastic Virus Database"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 0)">Name ↓</th>
                    <th title="The date the model was made"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 1)">Date made ↓</th>
                    <th title="Start-stop of the model in the sequence"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 2)">Range ↓</th>
                    <th title="Length of the model"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 3)">Length ↓</th>
                    <th title="The average pLDDT score of the model"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 4)">Average pLDDT ↓</th>
                    <th title="The H-score of the model"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 5)">H-score ↓</th>
                    <th title="The sequence identity of the model to the target protein"
                        onclick="sortTableAndFeatures('bfvd_table', sequence, ft4, BFVD_rowFeatureMap, '#div5', 6)">Sequence Identity ↓</th>
                </tr>
            </thead>
            <tbody>
                <tr v-for="bfvdmodel in bfvd_models" :data-feature-name="bfvdmodel.model_id">
                    <td><a v-bind:href="bfvdmodel.model_url" target="_blank">{{ bfvdmodel.model_id }}</a></td>
                    <td>{{ bfvdmodel.date_made }}</td>
                    <td>{{ bfvdmodel.range }}</td>
                    <td>{{ bfvdmodel.length }}</td>
                    <td>{{ bfvdmodel.avg_plddt | decimalPlaces }}</td>
                    <td>{{ bfvdmodel.h_score }}</td>
                    <td>{{ bfvdmodel.seq_ident }}</td>
                </tr>
            </tbody>
        </table>
    </div>
    `
});

Vue.component('esm-table', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            esm_models: this.$root.esm_models,
        }
    },
    mounted() {
        let show = Object.keys(this.$root.esm_models).length !== 0;
        toggleDisplay("esm-title", show);
        toggleDisplay("div6", show);
      },
    template: `
    <div id="esm-table-container"> 
            <table id="esm_table">
                <thead>
                    <tr>
                        <th title="Name of the model in the ESMfold Atlas database"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 0)">Name ↓</th>
                        <th title="The date the model was made"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 1)">Date made ↓</th>
                        <th title="Start-stop of the model in the sequence"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 2)">Range ↓</th>
                        <th title="Length of the model"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 3)">Length ↓</th>
                        <th title="The average pLDDT score of the model"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 4)">Average pLDDT ↓</th>
                        <th title="The H-score of the model"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 5)">H-score ↓</th>
                        <th title="The sequence identity of the model to the target protein"
                            onclick="sortTableAndFeatures('esm_table', sequence, ft5, ESM_rowFeatureMap, '#div6', 6)">Sequence Identity ↓</th>
                    </tr>
                </thead>
                <tbody>
                    <tr v-for="esmmodel in esm_models" :data-feature-name="esmmodel.model_id">
                        <td><a v-bind:href="esmmodel.model_url" target="_blank">{{ esmmodel.model_id }}</a></td>
                        <td>{{ esmmodel.date_made }}</td>
                        <td>{{ esmmodel.range }}</td>
                        <td>{{ esmmodel.length }}</td>
                        <td>{{ esmmodel.avg_plddt | decimalPlaces }}</td>
                        <td>{{ esmmodel.h_score }}</td>
                        <td>{{ esmmodel.seq_ident }}</td>
                    </tr>
                </tbody>
            </table>
        </div>
    `
});

Vue.component('hklinfo-feature-viewer', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            ss_pred: this.$root.ss_pred,
            classification: this.$root.classification,
            features1: [],
        }
    },
    template: `
      <div id="hklininfo-feature-viewer-container" class="content"></div>
    `,
    methods: {
        addFeature(feature) {
            this.features1.push(feature);
            this.ft1.addFeature(feature);
        },
        generateSecondaryStructureFeatures() {
            if (this.ss_pred) {
                let feature = {
                    data: this.ss_pred.regions.map(region => ({
                        x: region.start,
                        y: region.end,
                        color: region.colour,
                        description: region.text,
                    })),
                    name: "SS classification",
                    fontsize: "12",
                    className: "test1",
                    type: "multipleRect",
                    filter: "type1",
                };
                this.addFeature(feature);
            }
        },
        generateClassificationFeatures() {
            if (this.classification) {
                let feature = {
                    data: this.classification.regions.map(region => ({
                        x: region.start,
                        y: region.end,
                        color: region.colour,
                        description: region.text,
                    })),
                    name: "Classification",
                    fontsize: "12",
                    className: "test1",
                    type: "multipleRect",
                    filter: "type1",
                };
                this.addFeature(feature);
            }
        },
    },
    mounted() {
        var options = {
            showAxis: true,
            showSequence: true,
            brushActive: true,
            toolbar:true,
            bubbleHelp:true,
            zoomMax:10,
        };

        let show = this.$root.ss_pred || this.$root.classification ? true : false;
        if (show) {
            this.ft1 = new FeatureViewer.createFeature(sequence,"#div2", options);
            this.generateSecondaryStructureFeatures();
            this.generateClassificationFeatures();
        } else {
            toggleDisplay("seqinfo-title", false);
        }

    collapseDiv("collapsible2");
    }
});

Vue.component('pdb-feature-viewer', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            homologs: this.$root.homologs,
            pdb_features: [],
        }
    },
    template: `
      <div id="pdb-feature-viewer-container" class="content"></div>
    `,
    methods: {
        addFeature(feature) {
            this.pdb_features.push(feature);
            ft2.addFeature(feature);
            var rowId = feature.filter;
            window.PDB_rowFeatureMap[rowId] = feature;
        },
        generatePDBFeatures() {
            this.homologs.map(pdbs => {
                let feature = {
                    data: [{
                        x: pdbs._pfam_json.regions[0].start,
                        y: pdbs._pfam_json.regions[0].end,
                        color: pdbs._pfam_json.regions[0].colour,
                        description: pdbs._pfam_json.regions[0].text,
                    }],
                    name: pdbs._pfam_json.regions[0].text,
                    className: "class-" + pdbs._pfam_json.regions[0].text,
                    type: "multipleRect",
                    filter: "PDB-" + pdbs._pfam_json.regions[0].text,
                };
                this.addFeature(feature);
            });
        }

    },
    mounted() {
        window.PDB_rowFeatureMap = {};
        var options = {
            showAxis: true,
            showSequence: true,
            brushActive: true,
            toolbar:true,
            bubbleHelp:true,
            zoomMax:10,
        };
        window.ft2 = new FeatureViewer.createFeature(sequence,"#div3", options);
        this.generatePDBFeatures();

    selectTableAndFeatures(ft2);
    collapseDiv("collapsible3");
    }
});

Vue.component('af-feature-viewer', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            af_models: this.$root.af_models,
            af_features: [],
        }
    },
    template: `
      <div id="af-feature-viewer-container" class="content"></div>
    `,
    methods: {
        addFeature(feature) {
            this.af_features.push(feature);
            ft3.addFeature(feature);
            var rowId = feature.filter;
            window.AF_rowFeatureMap[rowId] = feature;
        },
        generateAFFeatures() {
            const colors = {
                'v_low': '#FF7D45',
                'low': '#FFDB13',
                'confident': '#65CBF3',
                'v_high': '#0053D6'
            };

            const descriptions = {
                'v_low': 'Very Low Confidence (pLDDT < 50)',
                'low': 'Low Confidence (70 > pLDDT > 50)',
                'confident': 'Confident (90 > pLDDT > 70)',
                'v_high': 'Very High Confidence (pLDDT > 90)'
            };
        
            if (this.af_models) {
                this.af_models.forEach(model => {
                    const modelName = model.name;
                    let data = [];
                    for (const [confidence, regions] of Object.entries(model.plddt_regions)) {
                        regions.forEach(([start, end]) => {
                            data.push({
                                x: start,
                                y: end,
                                color: colors[confidence],
                                description: descriptions[confidence],
                            });
                        });
                    }
                    let feature = {
                        data: data,
                        name: modelName,
                        className: modelName,
                        type: "multipleRect",
                        filter: modelName,
                    };
                    this.addFeature(feature);
                });
            }
        }

    },
    mounted() {
        window.AF_rowFeatureMap = {};
        var options = {
            showAxis: true,
            showSequence: true,
            brushActive: true,
            toolbar:true,
            bubbleHelp:true,
            zoomMax:10,
        };
        window.ft3 = new FeatureViewer.createFeature(sequence,"#div4", options);
        this.generateAFFeatures();

    selectTableAndFeatures(ft3);
    collapseDiv("collapsible4");
    }
});

Vue.component('bfvd-feature-viewer', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            bfvd_models: this.$root.bfvd_models,
            bfvd_features: [],
        }
    },
    template: `
      <div id="bfvd-feature-viewer-container" class="content"></div>
    `,
    methods: {
        addFeature(feature) {
            this.bfvd_features.push(feature);
            ft4.addFeature(feature);
            var rowId = feature.filter;
            window.BFVD_rowFeatureMap[rowId] = feature;
        },
        generateBFVDFeatures() {
            const colors = {
                'v_low': '#FF7D45',
                'low': '#FFDB13',
                'confident': '#65CBF3',
                'v_high': '#0053D6'
            };

            const descriptions = {
                'v_low': 'Very Low Confidence (pLDDT < 50)',
                'low': 'Low Confidence (70 > pLDDT > 50)',
                'confident': 'Confident (90 > pLDDT > 70)',
                'v_high': 'Very High Confidence (pLDDT > 90)'
            };
        
            if (this.bfvd_models) {
                this.bfvd_models.forEach(model => {
                    const modelName = model.name;
                    let data = [];
                    for (const [confidence, regions] of Object.entries(model.plddt_regions)) {
                        regions.forEach(([start, end]) => {
                            data.push({
                                x: start,
                                y: end,
                                color: colors[confidence],
                                description: descriptions[confidence],
                            });
                        });
                    }
                    let feature = {
                        data: data,
                        name: modelName,
                        className: modelName,
                        type: "multipleRect",
                        filter: modelName,
                    };
                    this.addFeature(feature);
                });
            }
        }

    },
    mounted() {
        window.BFVD_rowFeatureMap = {};
        var options = {
            showAxis: true,
            showSequence: true,
            brushActive: true,
            toolbar:true,
            bubbleHelp:true,
            zoomMax:10,
        };
        window.ft4 = new FeatureViewer.createFeature(sequence,"#div5", options);
        this.generateBFVDFeatures();

    selectTableAndFeatures(ft4);
    collapseDiv("collapsible5");
    }
});

Vue.component('esm-feature-viewer', {
    data: function () {
        if (typeof window.sequence === 'undefined') {
            window.sequence = this.$root.sequence;
        }
        return {
            esm_models: this.$root.esm_models,
            esm_features: [],
        }
    },
    template: `
      <div id="esm-feature-viewer-container" class="content"></div>
    `,
    methods: {
        addFeature(feature) {
            this.esm_features.push(feature);
            ft5.addFeature(feature);
            var rowId = feature.filter;
            window.ESM_rowFeatureMap[rowId] = feature;
        },
        generateESMFeatures() {
            const colors = {
                'v_low': '#FF7D45',
                'low': '#FFDB13',
                'confident': '#65CBF3',
                'v_high': '#0053D6'
            };

            const descriptions = {
                'v_low': 'Very Low Confidence (pLDDT < 50)',
                'low': 'Low Confidence (70 > pLDDT > 50)',
                'confident': 'Confident (90 > pLDDT > 70)',
                'v_high': 'Very High Confidence (pLDDT > 90)'
            };
        
            if (this.esm_models) {
                this.esm_models.forEach(model => {
                    const modelName = model.name;
                    let data = [];
                    for (const [confidence, regions] of Object.entries(model.plddt_regions)) {
                        regions.forEach(([start, end]) => {
                            data.push({
                                x: start,
                                y: end,
                                color: colors[confidence],
                                description: descriptions[confidence],
                            });
                        });
                    }
                    let feature = {
                        data: data,
                        name: modelName,
                        className: modelName,
                        type: "multipleRect",
                        filter: modelName,
                    };
                    this.addFeature(feature);
                });
            }
        }

    },
    mounted() {
        window.ESM_rowFeatureMap = {};
        var options = {
            showAxis: true,
            showSequence: true,
            brushActive: true,
            toolbar:true,
            bubbleHelp:true,
            zoomMax:10,
        };
        window.ft5 = new FeatureViewer.createFeature(sequence,"#div6", options);
        this.generateESMFeatures();

    selectTableAndFeatures(ft5);
    collapseDiv("collapsible6");
    }
});


new Vue({
    el: '#app',
    data: {
        homologs: mrparse_data.pfam.homologs,
        ss_pred: mrparse_data.pfam.ss_pred,
        classification: mrparse_data.pfam.classification,
        hklinfo: mrparse_data.hkl_info,
        af_models: mrparse_data.pfam.af_models,
        bfvd_models: mrparse_data.pfam.bfvd_models,
        esm_models: mrparse_data.pfam.esm_models,
        sequence: mrparse_data.pfam.sequence,
    },
})
