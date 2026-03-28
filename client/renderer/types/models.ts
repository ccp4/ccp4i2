/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
export class Project {
  constructor(
    public tags: number[] | ProjectTag[],
    public exports: number[],
    public imports: number[],
    public jobs: number[],
    public id: number,
    public uuid: string,
    public name: string,
    public description: string,
    public directory: string,
    public creation_time: string,
    public creation_user: string,
    public creation_host: string,
    public last_access: string,
    public last_job_number: number,
    public follow_from_job: number,
    public i1_project_name: string,
    public i1_project_directory: string
  ) {}
}
export class ProjectTag {
  constructor(
    public children: number[],
    public id: number,
    public parent: number,
    public text: string,
    public projects: number[]
  ) {}
}
export type ProjectExportStatus = "pending" | "running" | "completed" | "failed";

export class ProjectExport {
  constructor(
    public id: number,
    public project: Project | number,
    public time: string,
    public status: ProjectExportStatus,
    public file_exists?: boolean
  ) {}
}
export class ProjectImport {
  constructor(
    public id: number,
    public project: number,
    public time: string
  ) {}
}
export class Job {
  constructor(
    //public children: number[],
    public serverjob: number,
    public float_values: number[],
    public char_values: number[],
    public files: number[],
    public file_uses: number[],
    public xdatas: number[],
    public id: number,
    public uuid: string,
    public project: number,
    public parent: number,
    public number: string,
    public title: string,
    public status: number,
    public evaluation: number,
    public comments: string,
    public creation_time: string,
    public finish_time: string,
    public task_name: string,
    public process_id: number,
    public directory?: string
  ) {}
}
export class ServerJob {
  constructor(
    public job: number,
    public server_process_id: number,
    public machine: string,
    public username: string,
    public mechanism: string,
    public remote_path: string,
    public custom_code_file: string,
    public validate: string,
    public key_file_name: string,
    public server_group: string
  ) {}
}
export class JobValueKey {
  constructor(
    public name: string,
    public description: string
  ) {}
}
export class JobFloatValue {
  constructor(
    public id: number,
    public job: number,
    public key: number,
    public value: number
  ) {}
}
export class JobCharValue {
  constructor(
    public id: number,
    public job: number,
    public key: number,
    public value: string
  ) {}
}
export class FileType {
  constructor(
    public files: number[],
    public name: string,
    public description: string
  ) {}
}
export class File {
  constructor(
    public exports: number[],
    public fileimport: number,
    public file_uses: number[],
    public id: number,
    public uuid: string,
    public name: string,
    public directory: number,
    public type: string,
    public sub_type: number,
    public content: number,
    public annotation: string,
    public job: number,
    public job_param_name: string,
    public path?: string,
  ) {}
}
export class FileExport {
  constructor(
    public id: number,
    public file: number,
    public time: string,
    public name: string
  ) {}
}
export class FileImport {
  constructor(
    public file: number,
    public time: string,
    public name: string,
    public checksum: string,
    public last_modified: string
  ) {}
}
export class FileUse {
  constructor(
    public id: number,
    public file: number,
    public job: number,
    public role: number,
    public job_param_name: string
  ) {}
}
export class XData {
  constructor(
    public id: number,
    public data_class: string,
    public xml: string,
    public job: number
  ) {}
}

/** Job node with embedded files and KPIs from job_tree endpoint.
 * Extends Job but overrides files from number[] (IDs) to full File objects.
 */
export interface JobTreeNode extends Omit<Job, 'files'> {
  files: File[];
  kpis: {
    float_values: Record<string, number>;
    char_values: Record<string, string>;
  };
  children: JobTreeNode[];
}

/** Response from projects/{id}/job_tree endpoint */
export interface JobTreeResponse {
  job_tree: JobTreeNode[];
  total_jobs: number;
  total_files: number;
}

export const nullFile = {
  exports: [],
  fileimport: -1,
  file_uses: [],
  id: -1,
  uuid: "",
  name: "",
  directory: -1,
  type: "",
  sub_type: -1,
  content: -1,
  annotation: "",
  job: -1,
  job_param_name: "",
};

/**
 * Standard Django REST Framework paginated response.
 */
export interface PaginatedResponse<T> {
  count: number;
  next: string | null;
  previous: string | null;
  results: T[];
}
