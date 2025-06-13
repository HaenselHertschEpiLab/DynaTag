# Galaxy workflow: DynaTag_Response_23_TP53_DynaTag_WT.ChIPseq_R248Q.ChIP.seq_in_TP53_Peaks.WT.cond.exp._Peaks.R248Q_Peaks.DynaTag
```Galaxy
{
    "a_galaxy_workflow": "true",
    "annotation": "",
    "comments": [],
    "format-version": "0.1",
    "name": "DynaTag_Response_23_TP53_DynaTag_WT.ChIPseq_R248Q.ChIP.seq_in_TP53_Peaks.WT.cond.exp._Peaks.R248Q_Peaks.DynaTag",
    "report": {
        "markdown": "\n# Workflow Execution Report\n\n## Workflow Inputs\n```galaxy\ninvocation_inputs()\n```\n\n## Workflow Outputs\n```galaxy\ninvocation_outputs()\n```\n\n## Workflow\n```galaxy\nworkflow_display()\n```\n"
    },
    "steps": {
        "0": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 0,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "bigwig_Pascal_PDX_merged"
                }
            ],
            "label": "bigwig_Pascal_PDX_merged",
            "name": "Input dataset collection",
            "outputs": [],
            "position": {
                "left": 0,
                "top": 220
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bigwig\"], \"tag\": null, \"collection_type\": \"list\"}",
            "tool_version": null,
            "type": "data_collection_input",
            "uuid": "98012971-ba5c-4490-8888-bc997aad1926",
            "when": null,
            "workflow_outputs": []
        },
        "1": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 1,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "SRR12884191_treated_WT_p53_ChIPseq.bigwig"
                }
            ],
            "label": "SRR12884191_treated_WT_p53_ChIPseq.bigwig",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 410,
                "top": 550
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bigwig\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "7c10990b-f641-485a-8f12-84529cb05762",
            "when": null,
            "workflow_outputs": []
        },
        "2": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 2,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "SRR12884191_treated_WT_p53_ChIPseq.narrowpeaks"
                }
            ],
            "label": "SRR12884191_treated_WT_p53_ChIPseq.narrowpeaks",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 710,
                "top": 0
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bed\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "630ab475-a568-4f7f-b9eb-f3bf1eef1367",
            "when": null,
            "workflow_outputs": []
        },
        "3": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 3,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "SRR27596890_nontreated_R248Q_p53_ChIPseq.narrowpeaks"
                }
            ],
            "label": "SRR27596890_nontreated_R248Q_p53_ChIPseq.narrowpeaks",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 710,
                "top": 110
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bed\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "86e92aac-2446-49c2-b9cc-ea74f75098c1",
            "when": null,
            "workflow_outputs": []
        },
        "4": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 4,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "P53_all_peaks.over59nt.sorted.bed"
                }
            ],
            "label": "P53_all_peaks.over59nt.sorted.bed",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 720,
                "top": 250
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bed\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "1b18d16d-97f9-4008-ae98-b986b2489a19",
            "when": null,
            "workflow_outputs": []
        },
        "5": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 5,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "SRR12884190_nontreated_WT_p53_ChIPseq.bigwig"
                }
            ],
            "label": "SRR12884190_nontreated_WT_p53_ChIPseq.bigwig",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 410,
                "top": 680
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bigwig\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "5dd6075a-dc23-47b8-90d8-ccb676744f8a",
            "when": null,
            "workflow_outputs": []
        },
        "6": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 6,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "SRR27596890_nontreated_R248Q_p53_ChIPseq.bigwig"
                }
            ],
            "label": "SRR27596890_nontreated_R248Q_p53_ChIPseq.bigwig",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 410,
                "top": 810
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bigwig\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "5dd14344-736d-4362-89f6-31f1c290a906",
            "when": null,
            "workflow_outputs": []
        },
        "7": {
            "annotation": "",
            "content_id": null,
            "errors": null,
            "id": 7,
            "input_connections": {},
            "inputs": [
                {
                    "description": "",
                    "name": "SRR27596891_Input_R248Q_p53_ChIPseq.bigwig"
                }
            ],
            "label": "SRR27596891_Input_R248Q_p53_ChIPseq.bigwig",
            "name": "Input dataset",
            "outputs": [],
            "position": {
                "left": 410,
                "top": 930
            },
            "tool_id": null,
            "tool_state": "{\"optional\": false, \"format\": [\"bigwig\"], \"tag\": null}",
            "tool_version": null,
            "type": "data_input",
            "uuid": "2a8dd971-a409-4431-9e60-0d8f2d9e5771",
            "when": null,
            "workflow_outputs": []
        },
        "8": {
            "annotation": "",
            "content_id": "__EXTRACT_DATASET__",
            "errors": null,
            "id": 8,
            "input_connections": {
                "input": {
                    "id": 0,
                    "output_name": "output"
                }
            },
            "inputs": [],
            "label": "P53.CTRL.PDXS02730.cpm.bw",
            "name": "Extract dataset",
            "outputs": [
                {
                    "name": "output",
                    "type": "data"
                }
            ],
            "position": {
                "left": 400,
                "top": 200
            },
            "post_job_actions": {
                "ChangeDatatypeActionoutput": {
                    "action_arguments": {
                        "newtype": "bigwig"
                    },
                    "action_type": "ChangeDatatypeAction",
                    "output_name": "output"
                }
            },
            "tool_id": "__EXTRACT_DATASET__",
            "tool_state": "{\"input\": {\"__class__\": \"ConnectedValue\"}, \"which\": {\"which_dataset\": \"by_identifier\", \"__current_case__\": 1, \"identifier\": \"P53.CTRL.PDXS02730.cpm.bw\"}, \"__page__\": 0, \"__rerun_remap_job_id__\": null}",
            "tool_version": "1.0.1",
            "type": "tool",
            "uuid": "23a52da9-7ff4-47b8-9949-5b52444adfed",
            "when": null,
            "workflow_outputs": []
        },
        "9": {
            "annotation": "",
            "content_id": "__EXTRACT_DATASET__",
            "errors": null,
            "id": 9,
            "input_connections": {
                "input": {
                    "id": 0,
                    "output_name": "output"
                }
            },
            "inputs": [],
            "label": "P53.CHEM.PDXS02730.cpm.bw",
            "name": "Extract dataset",
            "outputs": [
                {
                    "name": "output",
                    "type": "data"
                }
            ],
            "position": {
                "left": 410,
                "top": 380
            },
            "post_job_actions": {
                "ChangeDatatypeActionoutput": {
                    "action_arguments": {
                        "newtype": "bigwig"
                    },
                    "action_type": "ChangeDatatypeAction",
                    "output_name": "output"
                }
            },
            "tool_id": "__EXTRACT_DATASET__",
            "tool_state": "{\"input\": {\"__class__\": \"ConnectedValue\"}, \"which\": {\"which_dataset\": \"by_identifier\", \"__current_case__\": 1, \"identifier\": \"P53.CHEM.PDXS02730.cpm.bw\"}, \"__page__\": 0, \"__rerun_remap_job_id__\": null}",
            "tool_version": "1.0.1",
            "type": "tool",
            "uuid": "10b73b1c-e865-4fc4-aec5-dd8b74ae217d",
            "when": null,
            "workflow_outputs": []
        },
        "10": {
            "annotation": "",
            "content_id": "toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_compute_matrix/deeptools_compute_matrix/3.5.4+galaxy0",
            "errors": null,
            "id": 10,
            "input_connections": {
                "multibigwig_conditional|multibigwig_repeats_0|bigwigfiles": {
                    "id": 8,
                    "output_name": "output"
                },
                "multibigwig_conditional|multibigwig_repeats_1|bigwigfiles": {
                    "id": 9,
                    "output_name": "output"
                },
                "multibigwig_conditional|multibigwig_repeats_2|bigwigfiles": {
                    "id": 1,
                    "output_name": "output"
                },
                "multibigwig_conditional|multibigwig_repeats_3|bigwigfiles": {
                    "id": 5,
                    "output_name": "output"
                },
                "multibigwig_conditional|multibigwig_repeats_4|bigwigfiles": {
                    "id": 6,
                    "output_name": "output"
                },
                "multibigwig_conditional|multibigwig_repeats_5|bigwigfiles": {
                    "id": 7,
                    "output_name": "output"
                },
                "regionsFiles_0|regionsFile": {
                    "id": 2,
                    "output_name": "output"
                },
                "regionsFiles_1|regionsFile": {
                    "id": 3,
                    "output_name": "output"
                },
                "regionsFiles_2|regionsFile": {
                    "id": 4,
                    "output_name": "output"
                }
            },
            "inputs": [
                {
                    "description": "runtime parameter for tool computeMatrix",
                    "name": "advancedOpt"
                }
            ],
            "label": null,
            "name": "computeMatrix",
            "outputs": [
                {
                    "name": "outFileName",
                    "type": "deeptools_compute_matrix_archive"
                }
            ],
            "position": {
                "left": 1000,
                "top": 210
            },
            "post_job_actions": {},
            "tool_id": "toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_compute_matrix/deeptools_compute_matrix/3.5.4+galaxy0",
            "tool_shed_repository": {
                "changeset_revision": "a60c359ec43c",
                "name": "deeptools_compute_matrix",
                "owner": "bgruening",
                "tool_shed": "toolshed.g2.bx.psu.edu"
            },
            "tool_state": "{\"advancedOpt\": {\"showAdvancedOpt\": \"yes\", \"__current_case__\": 1, \"binSize\": \"50\", \"sortRegions\": \"keep\", \"sortUsing\": \"mean\", \"averageTypeBins\": \"mean\", \"missingDataAsZero\": true, \"skipZeros\": false, \"minThreshold\": null, \"maxThreshold\": null, \"scale\": null, \"metagene\": false, \"transcriptID\": \"\", \"exonID\": \"\", \"transcript_id_designator\": \"\", \"blackListFileName\": {\"__class__\": \"RuntimeValue\"}}, \"custom_sample_labels_conditional\": {\"custom_labels_select\": \"Yes\", \"__current_case__\": 1, \"labels\": \"DT_PDX_CTRL DT_PDX_CHEM WT_treated_exp WT_CTRL_exp R248Q_nontreated R248Q_Input\"}, \"mode\": {\"mode_select\": \"reference-point\", \"__current_case__\": 1, \"referencePoint\": \"center\", \"nanAfterEnd\": false, \"beforeRegionStartLength\": \"2500\", \"afterRegionStartLength\": \"2500\"}, \"multibigwig_conditional\": {\"orderMatters\": \"Yes\", \"__current_case__\": 1, \"multibigwig_repeats\": [{\"__index__\": 0, \"bigwigfiles\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 1, \"bigwigfiles\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 2, \"bigwigfiles\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 3, \"bigwigfiles\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 4, \"bigwigfiles\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 5, \"bigwigfiles\": {\"__class__\": \"RuntimeValue\"}}]}, \"output\": {\"showOutputSettings\": \"yes\", \"__current_case__\": 1, \"saveMatrix\": false, \"saveSortedRegions\": false}, \"regionsFiles\": [{\"__index__\": 0, \"regionsFile\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 1, \"regionsFile\": {\"__class__\": \"RuntimeValue\"}}, {\"__index__\": 2, \"regionsFile\": {\"__class__\": \"RuntimeValue\"}}], \"__page__\": 0, \"__rerun_remap_job_id__\": null}",
            "tool_version": "3.5.4+galaxy0",
            "type": "tool",
            "uuid": "ab364216-3f9d-4e4d-8770-09db08f3ddd0",
            "when": null,
            "workflow_outputs": []
        },
        "11": {
            "annotation": "",
            "content_id": "toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_plot_heatmap/deeptools_plot_heatmap/3.5.4+galaxy0",
            "errors": null,
            "id": 11,
            "input_connections": {
                "matrixFile": {
                    "id": 10,
                    "output_name": "outFileName"
                }
            },
            "inputs": [
                {
                    "description": "runtime parameter for tool plotHeatmap",
                    "name": "matrixFile"
                }
            ],
            "label": "plotHeatmap_TP53_DT_ChIPseq_hg38_2500bp_black_and_white",
            "name": "plotHeatmap",
            "outputs": [
                {
                    "name": "outFileName",
                    "type": "png"
                }
            ],
            "position": {
                "left": 1380,
                "top": 180
            },
            "post_job_actions": {
                "DeleteIntermediatesActionoutFileName": {
                    "action_arguments": {},
                    "action_type": "DeleteIntermediatesAction",
                    "output_name": "outFileName"
                },
                "RenameDatasetActionoutFileName": {
                    "action_arguments": {
                        "newname": "plotHeatmap_TP53_DT_ChIPseq_hg38_2500bp_black_and_white"
                    },
                    "action_type": "RenameDatasetAction",
                    "output_name": "outFileName"
                }
            },
            "tool_id": "toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_plot_heatmap/deeptools_plot_heatmap/3.5.4+galaxy0",
            "tool_shed_repository": {
                "changeset_revision": "e4c7985d4585",
                "name": "deeptools_plot_heatmap",
                "owner": "bgruening",
                "tool_shed": "toolshed.g2.bx.psu.edu"
            },
            "tool_state": "{\"advancedOpt\": {\"showAdvancedOpt\": \"yes\", \"__current_case__\": 1, \"sortRegions\": \"descend\", \"sortUsing\": \"mean\", \"sortUsingSamples\": null, \"linesAtTickMarks\": false, \"averageTypeSummaryPlot\": \"mean\", \"plotType\": \"lines\", \"missingDataColor\": \"black\", \"colorMapRepeat\": [{\"__index__\": 0, \"colorMap\": \"binary\"}], \"alpha\": \"1.0\", \"colorList\": \"\", \"zMin\": \"0\", \"zMax\": \"\", \"yMin\": null, \"yMax\": null, \"xAxisLabel\": \"Distance\", \"yAxisLabel\": \"Peaks_Genes\", \"heatmapWidth\": \"7.5\", \"heatmapHeight\": \"25.0\", \"whatToShow\": \"heatmap and colorbar\", \"startLabel\": \"Upstream\", \"endLabel\": \"Downstream\", \"referencePointLabel\": \"Peak centre\", \"samplesLabel\": \"DT_PDX_CTRL DT_PDX_CHEM WT_treated_exp WT_CTRL_exp R248Q_nontreated R248Q_Input\", \"regionsLabel\": \"WT_NSCLC R248Q_NSCLC PDX_R248Q_SCLC\", \"plotTitle\": \"TP53 DynaTag and ChIP-seq in Peaks (WT_NSCLC R248Q_NSCLC PDX_R248Q_SCLC) \", \"legendLocation\": \"best\", \"labelRotation\": \"0\", \"perGroup\": true, \"used_multiple_regions\": {\"used_multiple_regions_options\": \"yes\", \"__current_case__\": 1}, \"clusterUsingSamples\": null}, \"matrixFile\": {\"__class__\": \"RuntimeValue\"}, \"output\": {\"showOutputSettings\": \"yes\", \"__current_case__\": 1, \"outFileFormat\": \"pdf\", \"dpi\": \"200\", \"saveMatrix\": false, \"saveSortedRegions\": false}, \"__page__\": 0, \"__rerun_remap_job_id__\": null}",
            "tool_version": "3.5.4+galaxy0",
            "type": "tool",
            "uuid": "7b368288-8345-43a1-8e4b-04fdaf59d114",
            "when": null,
            "workflow_outputs": [
                {
                    "label": "plotHeatmap_TP53_DT_ChIPseq_hg38_2500bp_black_and_white",
                    "output_name": "outFileName",
                    "uuid": "fd3d9801-f16a-4d80-98ec-bc7afc1d07fd"
                }
            ]
        }
```
    },
    "tags": [],
    "uuid": "de02a7ba-b929-48f0-937f-2c4f20bd735f",
    "version": 1
}
